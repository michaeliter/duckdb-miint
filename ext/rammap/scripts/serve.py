#!/usr/bin/env python3
"""Dev server with COOP/COEP headers required for SharedArrayBuffer (WASM threads).

Usage: python3 scripts/serve.py [port]
  Default port: 8080
  Open: http://localhost:8080/web/
"""

import sys
from http.server import HTTPServer, SimpleHTTPRequestHandler

class COEPHandler(SimpleHTTPRequestHandler):
    def end_headers(self):
        self.send_header('Cross-Origin-Opener-Policy', 'same-origin')
        self.send_header('Cross-Origin-Embedder-Policy', 'require-corp')
        # Disable browser caching during dev: rebuilds happen often and a
        # stale worker.js paired with a fresh index.html silently breaks the
        # demo.
        self.send_header('Cache-Control', 'no-store, max-age=0')
        super().end_headers()

port = int(sys.argv[1]) if len(sys.argv) > 1 else 8080
print(f'Serving at http://localhost:{port}/web/')
print('(COOP/COEP headers enabled for SharedArrayBuffer)')
HTTPServer(('0.0.0.0', port), COEPHandler).serve_forever()
