# Dev container

"Reopen in Container" works without any setup: it uses the shared
`devcontainer.json`. A developer who needs personal settings adds a local
config next to it, which git ignores.

| Path | Tracked | Contents |
|---|---|---|
| `devcontainer.json` | yes | Shared config. Builds `Dockerfile` directly. |
| `Dockerfile` | yes | Shared image: Python 3.12, uv, the C++ toolchain. |
| `post-create.sh` | yes | Shared first-run steps: `uv sync --all-extras`, the pre-commit hooks, a smoke test. |
| `local.devcontainer.json.example`, `local.Dockerfile.example` | yes | Templates of the local config. |
| `local/devcontainer.json`, `local/Dockerfile` | no | The developer's own config: mounts, resources, extra tools. |

## Local config

```bash
make devcontainer    # copies the templates to .devcontainer/local/; keeps existing files
```

Edit `local/`, then pick "quicksynth (local)" when VS Code asks which config to
open. Before each build, the local config builds the shared image on the host
as `quicksynth-dev-base` (`initializeCommand`). `local/Dockerfile` then starts
`FROM quicksynth-dev-base`. That works with Docker's default builder, which
keeps images in the local image store. A buildx builder with the
`docker-container` driver does not see them. A change to the templates does
not reach existing `local/` copies; compare them after pulling.

## Python environment

Each checkout and worktree gets its own environment, `.venv-container/` in its
root (`UV_PROJECT_ENVIRONMENT=.venv-container`, relative to the project root).
A host `.venv` and other checkouts in the same container never share it. VS
Code's interpreter and `PATH` point to the opened checkout's environment.
