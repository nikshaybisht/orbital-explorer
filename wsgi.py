"""WSGI entrypoint for production servers (gunicorn, etc.).

    gunicorn wsgi:server --bind 0.0.0.0:$PORT --workers 1 --threads 4
"""

from __future__ import annotations

import os

from orbital_explorer.app import build_app

_dash_app = build_app()
server = _dash_app.server

# Dash reads this when using its own server; gunicorn uses `server` above.
application = server

if __name__ == "__main__":
    port = int(os.environ.get("PORT", "8050"))
    host = os.environ.get("HOST", "0.0.0.0")
    _dash_app.run(debug=False, host=host, port=port)
