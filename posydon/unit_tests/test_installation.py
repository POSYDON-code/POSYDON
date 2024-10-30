import subprocess

def test_installation():
    """Test that POSYDON installs successfully."""
    result = subprocess.run(["pip", "install", "."], capture_output=True, text=True)
    assert result.returncode == 0, f"Installation failed: {result.stderr}"
