# Dev container

"Reopen in Container" works without any setup: it uses the shared
`devcontainer.json`. For personal settings, add a local config next to it,
which git ignores.

| Path | Tracked | Contents |
|---|---|---|
| `devcontainer.json` | yes | Shared config. Builds `Dockerfile` directly. |
| `Dockerfile` | yes | Shared image: Python 3.12, uv, the C++ toolchain. |
| `post-create.sh` | yes | Shared first-run steps: `uv sync --all-extras`, the pre-commit hooks, a smoke test. |
| `local.devcontainer.json.example`, `local.Dockerfile.example` | yes | Templates of the local config. |
| `local/devcontainer.json`, `local/Dockerfile` | no | Your own config: mounts, resources, extra tools. |

## Local config

```bash
make devcontainer    # copies the templates to .devcontainer/local/; keeps existing files
```

Edit the files in `local/`, then pick "quicksynth (local)" when VS Code asks
which config to open.

Before each build, the local config builds the shared image on the host as
`quicksynth-dev-base` (its `initializeCommand`), and `local/Dockerfile` starts
`FROM quicksynth-dev-base`. This works with Docker's default builder, which keeps
images in the local image store. A buildx builder with the `docker-container`
driver does not see them.

`make devcontainer` never overwrites your `local/` copies, so later changes to
the templates do not reach them. Compare them with the templates after pulling.

## Python environment

Each checkout and worktree gets its own environment, `.venv-container/` in its
root. The container sets `UV_PROJECT_ENVIRONMENT=.venv-container`, and uv
resolves the relative path against the project root. A host `.venv` and the
other checkouts in the container therefore never share it. An absolute path
would give them all one environment, and so one compiled `_core`. VS Code's
interpreter and `PATH`, in login shells too, point to the opened checkout's
environment.
