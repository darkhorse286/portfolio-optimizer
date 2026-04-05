import io
import builtins
import pytest
from scripts.generate_report import build_html_report


class _FakeWriter:
    def __init__(self):
        self._chunks = []
        self.closed = False

    def write(self, s):
        self._chunks.append(s)

    def getvalue(self):
        return "".join(self._chunks)

    def close(self):
        self.closed = True
    
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()
        return False


class _FakeOpen:
    def __init__(self, files):
        # files: dict path->bytes or str
        self._files = files

    def __call__(self, path, mode='r', *args, **kwargs):
        # Return a BytesIO for binary reads, FakeWriter for text writes
        data = self._files.get(path)
        if 'b' in mode:
            if data is None:
                raise FileNotFoundError(path)
            return io.BytesIO(data)
        else:
            # For write mode, create buffer to capture output
            if 'w' in mode:
                buf = _FakeWriter()
                # attach to dict so test can inspect
                self._files[path] = buf
                return buf
            # text read
            if data is None:
                raise FileNotFoundError(path)
            return io.StringIO(data.decode('utf-8') if isinstance(data, bytes) else str(data))


def test_build_html_report_missing_keys(monkeypatch):
    metrics = {'total_return': 0.05}  # many keys missing -> should show N/A

    # Prepare fake images as bytes
    equity_png = 'equity.png'
    drawdown_png = 'drawdown.png'
    fake_files = {
        equity_png: b"PNGDATA",
        drawdown_png: b"PNGDATA",
    }

    fake_open = _FakeOpen(fake_files)
    monkeypatch.setattr(builtins, 'open', fake_open)

    out_path = 'report.html'
    build_html_report(metrics, equity_png, drawdown_png, out_path, title='T')

    # Inspect written HTML
    written = fake_files[out_path].getvalue()
    assert 'Portfolio Report' in written or 'T' in written
    assert 'Total return' in written
    assert 'N/A' in written
    assert 'data:image/png;base64' in written
