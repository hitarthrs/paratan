"""Browser transport extension; uses Trame's Vue slot and VTK sync context."""
from pathlib import Path
import hashlib

static_path = Path(__file__).with_name('static')
# Hash both scripts so a dropzone change busts the browser cache too.
_asset_fingerprint = hashlib.sha256(
    (static_path / 'transport.js').read_bytes() + (static_path / 'dropzone.js').read_bytes()
).hexdigest()[:12]
asset_key = 'paratan-viewer-' + _asset_fingerprint
serve = {asset_key: str(static_path)}
scripts = [f'{asset_key}/transport.js', f'{asset_key}/dropzone.js']
vue_use = ['paratan_transport']
