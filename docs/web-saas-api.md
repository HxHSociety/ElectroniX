# Web SaaS API contract

The React app talks to a `Backend` (`frontend/src/api/`). On desktop that is
Tauri IPC; in the cloud it is `httpBackend`, which expects the server below.
**The server does not exist yet.** This document is the contract it must meet,
derived from what the frontend actually calls today.

Build the web app against a server with:

```bash
cd frontend
VITE_API_URL=https://api.example.com npm run build
```

## Transport

| Purpose | Request | Response |
|---|---|---|
| Run a command | `POST /rpc/{command}` JSON body = args | `200` JSON result, `204` no body, or `4xx/5xx {"error": "message"}` |
| Server push | `GET /events` (Server-Sent Events) | named events: `event: solver://done` + JSON `data:` |
| Fetch a file | `GET /assets?path=...` | file bytes (GLB models, etc.), `Range` supported |
| Upload a file | `PUT /uploads?name=...` raw body | `{"path": "<server path>"}` |

Authentication is a same-site **HttpOnly session cookie** (the client sends
`credentials: 'include'`). No tokens in JavaScript or URLs.

## Commands the frontend uses today

Verified against call sites in `frontend/src`:

| Command | Args | Returns |
|---|---|---|
| `pick_file` | (client-side: file chooser + `PUT /uploads`) | server path or `null` |
| `new_session` | - | `{ sessionId, ...dirs }` |
| `get_startup_session` | - | last session or `null` |
| `get_session_status` | `{ sessionId }` | `{ modelGlb, pcprep, thermalCsv, fatigueJson, creepPcres }` |
| `run_import` | `{ cvgPath, sessionId }` | GLB path |
| `run_generate_pc` | see `lib.rs` | `{ configPath, ... }` |
| `write_pcsim_file` | see `lib.rs` | `.pcsim` path |
| `run_solver` | `{ simId, pcsimPath, sessionId, simName }` | - (result arrives via `solver://done`) |
| `run_solver_auto` | see `lib.rs` | pcprep path |
| `read_json_file` | `{ path }` | file contents (string) |
| `read_trace_map` | `{ pcprepPath }` | trace nodes |
| `read_pcprep_all_nodes` | `{ pcprepPath, maxNodes }` | nodes |
| `read_pcres` | `{ path }` | creep records |

Events: `solver://done` (`{ simId, sessionId, fatiguePath?, pcprepPath?, error? }`)
and the progress event consumed in `SimulationBrowser.tsx`
(`{ step, line, done, error? }`).

Desktop-only, not implemented server-side: `check_for_update`
(the hook skips it unless the backend is `tauri`).

## Security: this is the part that must not be copied from the desktop app

The desktop commands take **absolute filesystem paths** from the UI
(`read_json_file { path }`, `read_pcres { path }`, `run_solver { pcsimPath }`).
That is fine for a local app and **a critical vulnerability on a server**: a
client could ask for `/etc/passwd`, another tenant's files, or cloud
credentials.

The server must instead:

1. Treat every `path` / `*Path` argument as an **opaque handle or a path
   relative to the caller's own session directory**, never as a host path.
2. Canonicalize and verify it resolves **inside** `sessions/<tenant>/<session>/`
   (reject `..`, symlinks that escape, absolute paths).
3. Authorize the session id against the authenticated user on every call.
4. Never pass client strings to a shell; sidecars are launched with an argument
   vector, and file names from uploads are sanitized.
5. Run the sidecars (`gltf_convertor`, `rpim_pc`, `rpim_solver`) in an
   unprivileged, network-less, resource-limited container per job.

## Suggested server shape

- `axum` API reusing the logic in `src-tauri/src/lib.rs`, ideally extracted into
  a shared crate so desktop and server cannot drift.
- Object storage (S3/GCS) under a per-tenant prefix; multipart/resumable
  uploads straight to storage for large CVG files (233 MB seen in testing).
- Queue + worker containers for `run_import`, `run_generate_pc`, `run_solver`.
  These are long jobs: return immediately, report progress via `/events`.
