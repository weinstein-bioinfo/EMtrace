# EMtrace on Windows

**Status:** experimental — not yet tested on a Windows machine.
**Last edited:** 2026-11-11

## Overview

EMtrace runs on Windows as a native application. Most of its external tools
have no Windows builds, so they run inside WSL (Windows Subsystem for Linux),
installed from Bioconda as on macOS and Linux. EMtrace finds them there
automatically; no configuration is needed for the standard setup below.

| Tool | Where it runs |
|:--|:--|
| EMtrace (GUI and CLI) | Windows |
| minimap2, samtools, bcftools, tabix | WSL |
| chopper, porechop_abi, pigz | WSL |
| DIAMOND, BLAST+ | WSL, or native Windows if installed there |
| Medaka (optional) | WSL, in its own environment |
<br>

For each tool EMtrace looks on the Windows `PATH` first and then inside WSL.
Paths are converted automatically (`C:\data\run1` becomes `/mnt/c/data/run1`).
Shell pipelines (read filtering, BAM sorting, variant calling) always run in
WSL.

The per-reference coverage plot is not drawn on Windows, because the Python
library it uses (*pysam*) has no Windows build. All other outputs are the
same.

## 1. Requirements

- Windows 10 version 2004 (build 19041) or later, or Windows 11
- 64-bit (x86_64) processor
- 16 GB RAM recommended; host filtering against the human genome needs more
  (see [Memory](#memory))
- Disk space for the databases, as on other platforms

## 2. Install EMtrace

1. Download `EMtrace-win-x86_64-{version}.zip` from the
   [GitHub releases](https://github.com/weinstein-bioinfo/EMtrace/releases).
2. Extract it, for example to `C:\Tools\EMtrace-x86_64\`.
3. Start `EMtrace.exe` in the extracted folder.

The app is not code-signed. If Windows shows "Windows protected your PC",
click **More info** → **Run anyway**.

<a id="enable-wsl"></a>
## 3. Enable WSL

Open **PowerShell as Administrator** and run:

```powershell
wsl --install
```

This installs WSL 2 with Ubuntu. Restart when prompted, then open **Ubuntu**
from the Start menu once to create your Linux user.

If WSL was already installed without a distribution:

```powershell
wsl --install -d Ubuntu
```

Check that it works:

```powershell
wsl -e echo "WSL OK"
```

## 4. Install the tools in WSL

Open **Ubuntu** from the Start menu and run:

```bash
curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b
~/miniforge3/bin/mamba create -y -n emtrace -c conda-forge -c bioconda \
    minimap2 samtools bcftools htslib chopper porechop_abi pigz diamond blast
```

The environment must be named `emtrace`: EMtrace looks for
`~/miniforge3/envs/emtrace/bin` (and the same folder under `mambaforge`,
`miniconda3`, `anaconda3` and `micromamba`). For any other location, see
[Settings](#settings).

Do not install the tools with `apt`. Ubuntu does not package chopper or
porechop_abi, and its samtools and bcftools are older than Bioconda's, so
variant calling may give different results from the macOS and Linux
installations.

Check from PowerShell:

```powershell
wsl -e bash -c "~/miniforge3/envs/emtrace/bin/minimap2 --version"
```

### Medaka (optional)

Medaka pins its own Python stack, so install it in a separate environment:

```bash
~/miniforge3/bin/mamba create -y -n medaka -c conda-forge -c bioconda medaka
```

Then set **assembly.medaka.path** in EMtrace to the Linux path of the
executable, for example `/home/<user>/miniforge3/envs/medaka/bin/medaka`.
A path starting with `/` is looked up inside WSL.

<a id="settings"></a>
## 5. Settings

Two optional settings are read from the environment, or else from
`%APPDATA%\EMtrace\user_config.yaml`:

| Environment variable | `user_config.yaml` key | Meaning |
|:--|:--|:--|
| `EMTRACE_WSL_DISTRO` | `wsl_distro` | WSL distribution with the tools (default: the WSL default distribution) |
| `EMTRACE_WSL_BIN` | `wsl_bin_dir` | Linux folder(s) holding the tools, separated by `:`, searched first |
<br>

Example `user_config.yaml` entries:

```yaml
wsl_distro: Ubuntu-24.04
wsl_bin_dir: /opt/bio/bin
```

Set `wsl_distro` if `wsl -l -v` shows another default distribution, for
example `docker-desktop`. Restart EMtrace after changing either setting.

## 6. Data location and speed

WSL 2 reads Windows drives (`C:\`, `/mnt/c`) more slowly than its own Linux
file system. This matters most for large files: host-filter indexes (`.mmi`)
and big FASTQ files.

EMtrace accepts paths inside WSL in their Windows form,
`\\wsl.localhost\Ubuntu\home\<user>\...`, and converts them to `/home/<user>/...`
for the tools. Keeping the database folder there may be faster; this has not
been benchmarked yet.

<a id="memory"></a>
### Memory

WSL 2 may use at most half of the computer's memory by default. To raise the
limit, create `C:\Users\<YourName>\.wslconfig`:

```ini
[wsl2]
memory=24GB
```

then run `wsl --shutdown` and start EMtrace again.

## 7. Troubleshooting

### "Required tool not found: minimap2"

- Check that WSL works: `wsl -e echo test`
- Check that the tool is installed:
  `wsl -e bash -c "ls ~/miniforge3/envs/emtrace/bin/minimap2"`
- If the environment is elsewhere, set `wsl_bin_dir` (see [Settings](#settings)).
- If you have several distributions, set `wsl_distro`.

### "needs WSL (Windows Subsystem for Linux), which was not found"

WSL or its Linux distribution is missing. Follow [Enable WSL](#enable-wsl).

### The first step takes a few seconds to start

WSL starts its Linux virtual machine on first use after login. Later tool
calls start quickly.

### Force Stop

Force Stop stops the running Linux tool through WSL. If a tool keeps running
afterwards, run `wsl --shutdown` in PowerShell.
