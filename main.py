"""
Entry point for the PyInstaller-bundled executable.
For development, run app.py directly instead.

When frozen, this file handles two modes:
  - Executor mode: called by app.py's subprocess with [exe, executor.py, pickle_path]
  - Server mode:   normal startup, opens browser and starts Flask
"""
import os
import sys

# Windowed builds (console=False) leave these as None, and multiprocessing workers
# write their tracebacks to sys.stderr, which would crash them.
if sys.stdout is None:
    sys.stdout = open(os.devnull, 'w')
if sys.stderr is None:
    sys.stderr = open(os.devnull, 'w')

import multiprocessing

# Must run before app/biaxial_analysis are imported: on Windows every Pool worker
# re-launches this exe, and freeze_support() turns this process into that worker.
multiprocessing.freeze_support()

if getattr(sys, 'frozen', False):
    # Detect executor mode: app.py spawns [sys.executable, _MEIPASS/executor.py, pickle_path]
    if len(sys.argv) >= 3 and os.path.basename(sys.argv[1]) == 'executor.py':
        import pickle
        from biaxial_analysis import run
        with open(sys.argv[2], 'rb') as f:
            obj = pickle.load(f)
        if len(sys.argv) > 3:
            run(obj, report_file=sys.argv[3])
        else:
            run(obj)
        sys.exit(0)

    # Server mode: chdir so app.py's relative paths (LICENSE, images/, tmp/) resolve correctly
    os.chdir(sys._MEIPASS)  # type: ignore[attr-defined]

import signal
import threading
import webbrowser
from app import app

PORT = 5444

app.jinja_env.globals['is_bundled'] = getattr(sys, 'frozen', False)


@app.route('/quit')
def quit_app():
    threading.Timer(0.5, lambda: os.kill(os.getpid(), signal.SIGTERM)).start()
    return '<p style="font-family:sans-serif;padding:20px;">Saka Beton has stopped. You can close this tab.</p>'


if __name__ == '__main__':
    if getattr(sys, 'frozen', False):
        threading.Timer(1.5, webbrowser.open, args=[f'http://127.0.0.1:{PORT}']).start()
    app.run(host='127.0.0.1', port=PORT, debug=False)
