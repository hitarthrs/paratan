"""Browser transport extension; uses Trame's Vue slot and VTK sync context."""
from pathlib import Path
import hashlib

static_path = Path(__file__).with_name('static')
asset_key = 'paratan-viewer-' + hashlib.sha256((static_path/'transport.js').read_bytes()).hexdigest()[:12]
serve = {asset_key: str(static_path)}
scripts = [f'{asset_key}/transport.js']
vue_use = ['paratan_transport']
