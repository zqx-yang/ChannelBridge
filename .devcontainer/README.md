# Developing ChannelBridge in a dev container

The container holds everything needed to develop and test the package, plus
PyMOL for the figure stages. **CAVER is not included** — see below.

## Running it

You need Docker and an editor that speaks the [Dev Containers
spec](https://containers.dev/). Three routes, in rough order of least effort:

1. **GitHub Codespaces** — needs nothing installed locally. Open the repo on
   GitHub → *Code* → *Codespaces* → *Create codespace*. Note that `~/.claude`
   survives stop/start but not a rebuild, which the volume mount in
   `devcontainer.json` handles.
2. **VS Code locally** — install Docker Desktop (or Docker Engine) and the
   [Dev Containers
   extension](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers),
   open the repo, and accept **Reopen in Container**. Or run *Dev Containers:
   Reopen in Container* from the Command Palette (`Ctrl/Cmd+Shift+P`).
3. **Dev Containers CLI** — `npm install -g @devcontainers/cli`, then
   `devcontainer up --workspace-folder .`.

On first build, `postCreateCommand` does an editable install with the dev
extras and runs the test suite, so a green build means the toolchain works.

Then open a terminal in the container and run `claude` to sign in.

## What is in the image

| | where it comes from | why |
|---|---|---|
| Python 3.11 + pip deps | `pip install -e '.[dev,figures,yaml]'` | the package and its test/lint tools |
| PyMOL (open source) | `apt install pymol` | the pose and figure stages **shell out** to it rather than importing it, so any 2.x/3.x build works |
| JRE (headless) | `apt install default-jre-headless` | CAVER is pure Java |
| Claude Code | the official dev container feature | |

`CHBRIDGE_PYMOL` is set to `/usr/bin/pymol`. The tool resolves PyMOL in the
order `CHBRIDGE_PYMOL` → `pymol` on `PATH` → the author's cluster path, so on
a host with its own PyMOL you only need to set that one variable.

## CAVER is not in the image, and cannot be

CAVER 3.0 is GPL-3 software from Masaryk University, distributed from
[caver.cz](https://caver.cz/) under a registration that we are not entitled to
redistribute around. Download it yourself (free for academic use) and put it
where the container expects:

```bash
# from the repository root, inside or outside the container
mkdir -p vendor
unzip ~/Downloads/caver_3.0.3.zip -d vendor/
mv vendor/caver_3.0.3/caver vendor/caver     # the INNER caver/ dir, the one
                                             # holding caver.jar and lib/
```

`$CAVER` is already pointed at `vendor/caver` by `devcontainer.json`. Add
`vendor/` to `.gitignore` — it is not ours to commit.

Any of Java 8, 11 or 21 runs it; 21 is what the image provides.

## Running the test suite

```bash
pytest -q          # the whole suite
ruff check src tests
```

Tests that need external data — the CAVER cross-build regression and the
iMARS verification — **skip themselves when the data is absent**, so a clean
checkout passes without CAVER, without the reference structures and without a
cluster. The PyMOL-dependent code is checked by its own selftests rather than
by pytest, because it must run inside PyMOL's interpreter:

```bash
pymol -cq src/channelbridge/render/run_pse.py -- --selftest
pymol -cq src/channelbridge/render/run_hotspots.py -- --selftest
```

Judge those by the printed `PASS`/`FAILED` marker, not the exit code: PyMOL
is prone to segfaulting during interpreter teardown after the test has
already finished. Run them with `PYTHONUNBUFFERED=1` or the marker can be
lost in the buffer when that happens.

## Caveats worth knowing

* A dev container is not a security boundary against a malicious repository.
  Don't mount host secrets (`~/.ssh`, cloud credential files) into it.
* The reference container in `anthropics/claude-code` adds an egress firewall
  (`init-firewall.sh`, plus `NET_ADMIN`/`NET_RAW` in `runArgs`). This config
  deliberately leaves that out — it is not required for Claude Code, and
  ChannelBridge's own network use is limited to PDB downloads. Add it if your
  threat model wants it.
* Ray-tracing in a container is CPU-only and no faster than on the host; the
  large figure renders are still minutes to tens of minutes.
