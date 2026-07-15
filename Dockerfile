# Orbital Explorer  -  production image
# Build:  docker build -t orbital-explorer .
# Run:    docker run --rm -p 8050:8050 orbital-explorer

FROM python:3.12-slim-bookworm

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PORT=8050 \
    HOST=0.0.0.0

# System libs used by scipy / scikit-image / matplotlib / RDKit wheels
RUN apt-get update && apt-get install -y --no-install-recommends \
        libxrender1 \
        libxext6 \
        libsm6 \
        libglib2.0-0 \
        libgomp1 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY requirements.txt pyproject.toml README.md LICENSE ./
COPY src ./src
COPY wsgi.py run_app.py ./

# Image serves the Dash app via gunicorn (not Streamlit).
# Drop streamlit from the install to keep the image smaller.
RUN pip install --upgrade pip \
    && pip install \
        "numpy>=1.26,<3" "scipy>=1.11,<2" "scikit-image>=0.22,<0.27" \
        "matplotlib>=3.8,<4" "rdkit>=2024.3.1" "plotly>=5.18,<7" \
        "dash>=2.17,<3" "requests>=2.31,<3" "gunicorn>=22,<27" \
    && pip install -e .

EXPOSE 8050

# Single worker: Dash keeps in-memory callback state; multi-worker needs sticky sessions.
CMD gunicorn wsgi:server \
    --bind 0.0.0.0:${PORT} \
    --workers 1 \
    --threads 4 \
    --timeout 120 \
    --access-logfile - \
    --error-logfile -
