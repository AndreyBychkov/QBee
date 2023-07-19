import subprocess


def build_html():
    subprocess.run("sphinx-build -b html docs/source docs/build")
