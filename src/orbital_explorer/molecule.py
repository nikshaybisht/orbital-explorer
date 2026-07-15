"""Turn whatever the user types into a classified RDKit molecule.

Accepts a SMILES string ("c1ccccc1"), a common name ("benzene"), or a name that
OPSIN understands ("buta-1,3-diene"). It then auto-classifies the species so the
rest of the tool knows what kind of MO treatment makes sense:

    atom        single atom            -> show atomic orbitals
    diatomic    two atoms              -> (diatomic MO diagram, future)
    pi-system   conjugated p-orbitals  -> Huckel MO treatment
    sigma-only  no conjugated system   -> no pi MO diagram to draw

Charge and radical electrons are detected so the pi-electron count (and hence
HOMO/LUMO filling) comes out right for cations, anions and radicals.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

from rdkit import Chem, rdBase
from rdkit.Chem import rdMolDescriptors

# Bare "annulene" is a class name; default it to the prototypical aromatic one.
_DEFAULT_ANNULENE = 18
_ANNULENE_RE = re.compile(r"^\[?(\d+)\]?[\s-]*annulene$")

# Tokens that are valid SMILES for something else but, in an MO tool, almost
# certainly mean the diatomic (CO = methanol as SMILES; NO = hydroxylamine).
_AMBIGUOUS_DIATOMIC = {"co", "no"}

# A few everyday names so the tool works offline even without OPSIN/Java.
_COMMON_NAMES: dict[str, str] = {
    "hydrogen": "[H][H]", "h2": "[H][H]",
    "oxygen": "O=O", "o2": "O=O",
    "nitrogen": "N#N", "n2": "N#N",
    "ethylene": "C=C", "ethene": "C=C",
    "acetylene": "C#C", "ethyne": "C#C",
    "butadiene": "C=CC=C", "1,3-butadiene": "C=CC=C",
    "benzene": "c1ccccc1",
    "naphthalene": "c1ccc2ccccc2c1",
    "pyridine": "c1ccncc1",
    "pyrrole": "c1cc[nH]c1",
    "furan": "c1ccoc1",
    "allyl": "[CH2]C=C", "allyl radical": "[CH2]C=C",
    "allyl cation": "[CH2+]C=C", "allyl anion": "[CH2-]C=C",
    "cyclobutadiene": "C1=CC=C1",
    "cyclooctatetraene": "C1=CC=CC=CC=C1", "cot": "C1=CC=CC=CC=C1",
    "cyclopentadienyl": "[CH-]1C=CC=C1",
    "tropylium": "[CH+]1C=CC=CC=C1",
    "azulene": "C1=CC=C2C=CC=C2C=C1", "anthracene": "c1ccc2cc3ccccc3cc2c1",
    "styrene": "C=Cc1ccccc1", "fulvene": "C=C1C=CC=C1",
    "hexatriene": "C=CC=CC=C", "1,3,5-hexatriene": "C=CC=CC=C",
    # Diatomics for the MO-diagram tab.
    "fluorine": "FF", "f2": "FF", "chlorine": "ClCl", "cl2": "ClCl",
    "hf": "F", "hydrogen fluoride": "F", "hcl": "Cl", "hydrogen chloride": "Cl",
    "co": "[C-]#[O+]", "carbon monoxide": "[C-]#[O+]",
    "no": "[N]=[O]", "nitric oxide": "[N]=[O]",
    "lih": "[Li][H]", "lithium hydride": "[Li][H]",
    # Noble-gas / exotic homonuclear diatomics: written as two separate atoms
    # because RDKit rejects e.g. a He-He bond. The diatomic engine only needs
    # the two element symbols, so a non-bonded pair is fine here.
    "he2": "[He].[He]", "ne2": "[Ne].[Ne]", "li2": "[Li].[Li]",
    "be2": "[Be].[Be]", "b2": "[B].[B]", "c2": "[C].[C]",
    "o2+": "[O+]=[O]", "o2-": "[O-]=[O]", "n2+": "[N+]#N",
    "water": "O", "methane": "C", "ethane": "CC",
}


@dataclass
class MoleculeInfo:
    query: str
    smiles: str | None = None
    formula: str | None = None
    charge: int = 0
    radical_electrons: int = 0
    n_heavy_atoms: int = 0
    category: str = "unknown"
    pi_atom_indices: list[int] = field(default_factory=list)
    note: str = ""
    mol: Chem.Mol | None = field(default=None, repr=False)

    @property
    def has_pi_system(self) -> bool:
        return self.category == "pi-system"


def _to_smiles(text: str) -> str | None:
    """Resolve a SMILES or name to a SMILES string.

    Order: SMILES as-is -> built-in name table (offline) -> OPSIN (offline, needs
    Java) -> PubChem REST (online). The last step is what makes most chemical
    names work; it is skipped gracefully when offline.
    """
    text = text.strip()
    if not text:
        return None
    key = text.lower()
    # "CO"/"NO" are valid SMILES (methanol / hydroxylamine) but in an MO tool the
    # user almost always means the diatomic; prefer the name for these.
    if key in _AMBIGUOUS_DIATOMIC:
        return _COMMON_NAMES[key]
    with rdBase.BlockLogs():
        if Chem.MolFromSmiles(text) is not None:
            return text
    annulene = _annulene_smiles(key)
    if annulene:
        return annulene
    if key in _COMMON_NAMES:
        return _COMMON_NAMES[key]
    try:
        import warnings

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # OPSIN warns loudly when Java is absent
            from py2opsin import py2opsin

            result = py2opsin(text)
        if result:
            return result.strip().splitlines()[0]
    except Exception:
        pass
    return _pubchem_smiles(text)


def _annulene_smiles(key: str) -> str | None:
    """[n]annulene -> a monocyclic fully-conjugated ring SMILES.

    "annulene" alone is a class name, so it defaults to [18]annulene (the
    prototypical aromatic annulene). "[n]annulene" / "n-annulene" build the
    n-membered ring. Odd or <4 ring sizes return None (no neutral closed shell).
    """
    if key == "annulene":
        n = _DEFAULT_ANNULENE
    else:
        match = _ANNULENE_RE.match(key)
        if not match:
            return None
        n = int(match.group(1))
    if n < 4 or n % 2 != 0:
        return None
    return "C1" + "=CC" * (n // 2 - 1) + "=C1"


def _pubchem_smiles(name: str) -> str | None:
    """Name -> SMILES via the PubChem REST API. Returns None on any failure
    (including being offline), so the offline paths above still work."""
    import urllib.parse

    try:
        import requests
    except ImportError:
        return None

    encoded = urllib.parse.quote(name)
    # PubChem renamed IsomericSMILES -> SMILES in 2025; try the new name first.
    for prop in ("SMILES", "IsomericSMILES"):
        url = (
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
            f"{encoded}/property/{prop}/JSON"
        )
        try:
            resp = requests.get(url, timeout=15)
            resp.raise_for_status()
            smiles = resp.json()["PropertyTable"]["Properties"][0].get(prop)
            if smiles:
                return smiles
        except Exception:
            continue
    return None


def _pi_atoms(mol: Chem.Mol) -> list[int]:
    """Indices of atoms that carry a p-orbital available for pi conjugation.

    Heuristic: aromatic or sp/sp2 atoms, plus any radical/charged centre that is
    bonded to such an atom (so the radical carbon of allyl is included).
    """
    sp = {Chem.HybridizationType.SP, Chem.HybridizationType.SP2}
    core = {
        a.GetIdx()
        for a in mol.GetAtoms()
        if a.GetIsAromatic() or a.GetHybridization() in sp
    }
    extra = set()
    for a in mol.GetAtoms():
        if a.GetIdx() in core:
            continue
        if a.GetNumRadicalElectrons() or a.GetFormalCharge():
            if any(nb.GetIdx() in core for nb in a.GetNeighbors()):
                extra.add(a.GetIdx())
    return sorted(core | extra)


def resolve(text: str) -> MoleculeInfo:
    """Resolve and classify a user-typed species into a MoleculeInfo."""
    info = MoleculeInfo(query=text)
    smiles = _to_smiles(text)
    if smiles is None:
        info.note = "Could not resolve that name. Try SMILES, check spelling, or check internet (names use PubChem)."
        return info

    with rdBase.BlockLogs():
        mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        info.note = f"RDKit could not parse SMILES {smiles!r}."
        return info

    info.mol = mol
    info.smiles = Chem.MolToSmiles(mol)
    info.formula = rdMolDescriptors.CalcMolFormula(mol)
    info.charge = Chem.GetFormalCharge(mol)
    info.radical_electrons = sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms())
    info.n_heavy_atoms = mol.GetNumHeavyAtoms()
    n_total = Chem.AddHs(mol).GetNumAtoms()

    if n_total == 1:
        info.category = "atom"
        info.note = "Single atom: use the Atomic orbitals tab."
        return info
    if n_total == 2:
        info.category = "diatomic"
        info.note = "Diatomic: use the MO diagram tab."
        return info

    pi = _pi_atoms(mol)
    if len(pi) >= 2:
        info.category = "pi-system"
        info.pi_atom_indices = pi
        info.note = _pi_note(mol, pi)
    else:
        info.category = "sigma-only"
        info.note = "No conjugated pi-system; there is no pi MO diagram to draw."
    return info


def _pi_components(mol: Chem.Mol, pi: list[int]) -> list[list[int]]:
    """Connected components of the pi-atom subgraph (atoms joined by a bond
    where both ends are pi atoms)."""
    members = set(pi)
    adj: dict[int, list[int]] = {i: [] for i in pi}
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if i in members and j in members:
            adj[i].append(j)
            adj[j].append(i)
    seen: set[int] = set()
    components = []
    for start in pi:
        if start in seen:
            continue
        stack, comp = [start], []
        while stack:
            n = stack.pop()
            if n in seen:
                continue
            seen.add(n)
            comp.append(n)
            stack.extend(adj[n])
        components.append(sorted(comp))
    return components


def _pi_note(mol: Chem.Mol, pi: list[int]) -> str:
    """An honest description: truly conjugated vs an isolated pi bond vs several
    separate pi units; warn about orthogonal pi systems (sp centres)."""
    comps = _pi_components(mol, pi)
    sizes = sorted((len(c) for c in comps), reverse=True)
    if len(comps) == 1 and len(pi) == 2:
        note = "Single pi bond (2 p-orbitals): bonding/antibonding pair."
    elif len(comps) == 1:
        note = f"Conjugated pi-system across {len(pi)} atoms; Huckel applies."
    else:
        note = f"{len(comps)} separate pi units (sizes {sizes}), not conjugated to each other."
    if any(mol.GetAtomWithIdx(i).GetHybridization() == Chem.HybridizationType.SP for i in pi):
        note += " Contains sp centres (allene/alkyne): the planar model shows only one of the orthogonal pi systems."
    return note


def is_monocyclic_pi(info: MoleculeInfo) -> bool:
    """True if the pi-system is exactly one ring (so the 4n+2 rule applies)."""
    if info.mol is None or not info.pi_atom_indices:
        return False
    rings = info.mol.GetRingInfo().AtomRings()
    return len(rings) == 1 and set(rings[0]) == set(info.pi_atom_indices)
