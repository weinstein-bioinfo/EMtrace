"""Database download and installation engine.

Port of ``install_db.pl``, adapted from virefset/utils/db_download.py.

Handles three database categories:

1. **Annotation databases** (``ann.databases`` in config.yaml) —
   minimap2 FASTA, BLASTN, Diamond targets.  All are tarballs
   hosted on CSC Allas (https://a3s.fi/lazypipe-db/).

2. **Host genome databases** (``host.databases`` in config.yaml) —
   minimap2-indexed FASTA tarballs on CSC Allas
   (https://a3s.fi/lazypipe-hostgen/).

3. **Taxonomy databases** — two sources:
   - NCBI: download taxdump.tar.gz → extract → build SQLite via
     ``NCBITaxonomy.build_database()``
   - SILVA: download taxmap .txt.gz → decompress in-place

Typical usage::

    from emtrace.utils.db_download import install_by_key
    from emtrace.core.config import EMtraceConfig

    cfg = EMtraceConfig.from_yaml("config.yaml")
    install_by_key(cfg, "minimap.nt.abv")
    install_by_key(cfg, "taxonomy")
    install_by_key(cfg, "Homo_sapiens")
    install_by_key(cfg, "hostdbs")   # all host genomes
"""

from __future__ import annotations

import datetime
import gzip
import logging
import os
import shutil
import tarfile
import urllib.request
from pathlib import Path
from typing import Callable, List, Optional, Tuple
from urllib.parse import urlparse

logger = logging.getLogger(__name__)

# ------------------------------------------------------------------
# Constants
# ------------------------------------------------------------------

_CHUNK_SIZE = 64 * 1024         # 64 KB read chunks
_TARBALL_SUFFIXES = (
    ".tar.gz", ".tgz", ".tar.bz2", ".tar.xz",
)

# Callback type: (bytes_done, bytes_total)
# bytes_total == -1 when Content-Length is unavailable
ProgressCallback = Optional[Callable[[int, int], None]]


# ------------------------------------------------------------------
# Exceptions
# ------------------------------------------------------------------

class DatabaseAlreadyInstalledError(Exception):
    """Database files already exist; use force=True to overwrite."""

    def __init__(self, key: str, path: Path) -> None:
        super().__init__(
            f"Database '{key}' already installed at {path}. "
            "Use force=True to overwrite."
        )
        self.key = key
        self.path = path


class DatabaseNotConfiguredError(Exception):
    """No URL or path is configured for this database."""

    def __init__(self, key: str) -> None:
        super().__init__(
            f"No URL configured for database '{key}'."
        )
        self.key = key


class DatabaseInstallError(Exception):
    """Download or extraction failed."""


# ------------------------------------------------------------------
# Low-level helpers
# ------------------------------------------------------------------

def _is_tarball(url: str) -> bool:
    path = urlparse(url).path
    return any(path.endswith(s) for s in _TARBALL_SUFFIXES)


def _is_plain_gz(url: str) -> bool:
    path = urlparse(url).path
    return path.endswith(".gz") and not _is_tarball(url)


def download_to_temp(
    url: str,
    dest_dir: Path,
    progress_cb: ProgressCallback = None,
) -> Path:
    """Stream *url* into a ``.tmp`` file under *dest_dir*.

    Supports ``http://``, ``https://``, and ``ftp://`` via
    ``urllib.request``.  The ``.tmp`` file is deleted on failure.

    Parameters
    ----------
    url : str
        Full URL to download.
    dest_dir : Path
        Directory for the temporary file.
    progress_cb : callable, optional
        Called as ``progress_cb(bytes_done, bytes_total)``
        where *bytes_total* is -1 when unknown.

    Returns
    -------
    Path
        Path to the downloaded ``.tmp`` file.
    """
    dest_dir.mkdir(parents=True, exist_ok=True)
    filename = Path(urlparse(url).path).name or "download"
    tmp_path = dest_dir / (filename + ".tmp")

    logger.info("Downloading %s", url)

    try:
        with urllib.request.urlopen(url) as resp:  # noqa: S310
            total = int(resp.headers.get("Content-Length", -1))
            done = 0
            with open(tmp_path, "wb") as fh:
                while True:
                    chunk = resp.read(_CHUNK_SIZE)
                    if not chunk:
                        break
                    fh.write(chunk)
                    done += len(chunk)
                    if progress_cb is not None:
                        progress_cb(done, total)
    except Exception:
        if tmp_path.exists():
            tmp_path.unlink()
        raise

    if total > 0 and done != total:
        tmp_path.unlink()
        raise DatabaseInstallError(
            f"Truncated download: got {done:,} of {total:,} "
            f"bytes from {url}"
        )

    size_mb = tmp_path.stat().st_size / 1e6
    logger.info(
        "Download complete: %s (%.1f MB)", tmp_path, size_mb,
    )
    return tmp_path


def extract_tarball(tmp_path: Path, target_dir: Path) -> None:
    """Extract *tmp_path* tarball into *target_dir*, then delete it."""
    target_dir.mkdir(parents=True, exist_ok=True)
    logger.info("Extracting %s → %s", tmp_path, target_dir)
    try:
        with tarfile.open(tmp_path) as tf:
            tf.extractall(target_dir)  # noqa: S202
    except EOFError:
        # Some archives lack a proper end-of-stream marker —
        # safe to ignore when all data was written.
        logger.debug(
            "Ignoring missing EOS marker in %s", tmp_path.name,
        )
    tmp_path.unlink()
    logger.info("Extraction complete: %s", target_dir)


def decompress_gz(tmp_path: Path, dest_path: Path) -> None:
    """Decompress plain ``.gz`` *tmp_path* to *dest_path*, delete tmp."""
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    logger.info("Decompressing %s → %s", tmp_path, dest_path)
    with gzip.open(tmp_path, "rb") as gz_in:
        with open(dest_path, "wb") as out:
            shutil.copyfileobj(gz_in, out)
    tmp_path.unlink()
    logger.info("Decompression complete: %s", dest_path)


def move_file(tmp_path: Path, dest_path: Path) -> None:
    """Atomically move *tmp_path* to *dest_path*."""
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path.replace(dest_path)
    logger.info("Installed: %s", dest_path)


# ------------------------------------------------------------------
# Installed-status helpers
# ------------------------------------------------------------------

def check_installed(db_path: Path) -> bool:
    """Return True if *db_path* or any file matching ``db_path*`` exists.

    Mirrors the Perl ``glob("$dbdir/$dbfile*")`` check.
    """
    if db_path.exists():
        return True
    return any(
        db_path.parent.glob(db_path.name + "*")
    )


def check_host_installed(db_path: Path) -> bool:
    """Return True if the host genome database is available on disk.

    The concatenated FASTA/gz is the database (minimap2 indexes on the fly),
    so the exact *db_path* file must exist.
    """
    return db_path.exists()


def list_installed_files(db_path: Path) -> List[Path]:
    """Return all files on disk matching ``db_path*``."""
    found = list(db_path.parent.glob(db_path.name + "*"))
    if db_path.exists() and db_path not in found:
        found.insert(0, db_path)
    return sorted(found)


def remove_database(db_path: Path) -> int:
    """Delete all files matching ``db_path*``.

    Returns the number of files deleted.
    """
    files = list_installed_files(db_path)
    for f in files:
        if f.is_dir():
            shutil.rmtree(f)
        else:
            f.unlink()
        logger.info("Removed: %s", f)
    return len(files)


# ------------------------------------------------------------------
# High-level install functions
# ------------------------------------------------------------------

def install_ann_database(
    key: str,
    url: str,
    db_path: Path,
    progress_cb: ProgressCallback = None,
    force: bool = False,
) -> None:
    """Download and install an annotation database.

    For tarball URLs the archive is extracted into the parent of
    *db_path*. For plain files the file is placed at *db_path*.

    Parameters
    ----------
    key : str
        Config key (used in error messages).
    url : str
        Download URL (from ``ann.databases[key].url``).
    db_path : Path
        Configured database path.
    progress_cb : callable, optional
        Progress callback ``(bytes_done, bytes_total)``.
    force : bool
        If False, raise if files already exist.
    """
    if not url or url in ("NA", "0"):
        raise DatabaseNotConfiguredError(key)

    if not force and check_installed(db_path):
        raise DatabaseAlreadyInstalledError(key, db_path)

    dest_dir = db_path.parent
    tmp = download_to_temp(url, dest_dir, progress_cb)

    if _is_tarball(url):
        extract_tarball(tmp, dest_dir)
    elif _is_plain_gz(url):
        decompress_gz(tmp, db_path)
    else:
        move_file(tmp, db_path)

    logger.info("Annotation database '%s' installed.", key)


def install_host_minimap2(
    key: str,
    urls: List[str],
    db_path: Path,
    progress_cb: ProgressCallback = None,
    force: bool = False,
) -> None:
    """Download and install a minimap2 host genome database.

    One or more FASTA/gz URLs are downloaded and concatenated into a single
    ``db_path`` file.  Concatenating gzip streams is valid and minimap2
    accepts the resulting multi-stream file without issues.
    """
    if not urls:
        raise DatabaseNotConfiguredError(key)
    if not force and check_host_installed(db_path):
        raise DatabaseAlreadyInstalledError(key, db_path)

    db_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_files: List[Path] = []
    try:
        for url in urls:
            tmp_files.append(download_to_temp(url, db_path.parent, progress_cb))
        # Concatenate all downloaded files into the target path
        with open(db_path, "wb") as out:
            for tmp in tmp_files:
                with open(tmp, "rb") as inp:
                    shutil.copyfileobj(inp, out)
    finally:
        for tmp in tmp_files:
            if tmp.exists():
                tmp.unlink()

    logger.info("Host genome (minimap2) '%s' installed: %s", key, db_path)


def _host_db_path(entry) -> Path:
    """Resolve a host ``DatabaseEntry`` to the on-disk path used at run time.

    Mirrors ``pipeline._resolve_host_databases`` exactly
    (``expandvars`` + ``expanduser`` of ``entry.db``) so that an index built
    here lands where ``BgFilterStep`` later looks for it.
    """
    return Path(os.path.expandvars(os.path.expanduser(entry.db)))


def _host_mmi_path(db_path: Path, k: int, w: int) -> Path:
    """Compute the ``.mmi`` index path for *db_path*.

    Matches ``BgFilterStep._align_minimap2``: a tuned ``k``/``w`` gets its own
    ``.k{k}w{w}.mmi`` suffix so it never clashes with the default ``.mmi``.
    """
    suffix = f".k{k}w{w}.mmi" if (k > 0 or w > 0) else ".mmi"
    return Path(str(db_path) + suffix)


def prebuild_host_index(
    cfg,
    key: str,
    *,
    force: bool = False,
    progress_cb: ProgressCallback = None,
) -> None:
    """Download (if needed) and build the minimap2 ``.mmi`` index for one host.

    Intended to be run serially **before** launching HPC array jobs, so that
    array tasks find a prebuilt index and never race to build it themselves
    (see ``Minimap2.host_filter`` / ``save_mmi_index``).

    Steps for the host database *key*:
      1. Ensure the FASTA is present (download + concatenate via
         :func:`install_host_minimap2` when missing, or always with *force*).
      2. Build the ``.mmi`` index at the path the run expects (derived from
         ``cfg.bgfilter.minimap.k/w``) unless it already exists (or *force*).
    """
    host_key = _find_host_key(cfg, key)
    if not host_key:
        raise DatabaseNotConfiguredError(key)
    entry = cfg.host_databases[host_key]

    db_path = _host_db_path(entry)

    # 1. FASTA present?
    if force or not check_host_installed(db_path):
        install_host_minimap2(
            host_key, entry.url_list, db_path,
            progress_cb=progress_cb, force=force,
        )
    else:
        logger.info("Host FASTA already present: %s", db_path)

    # 2. Build the index.
    mm_cfg = cfg.bgfilter.minimap
    mmi_path = _host_mmi_path(db_path, mm_cfg.k, mm_cfg.w)
    if mmi_path.exists() and not force:
        logger.info("Index already present, skipping: %s", mmi_path)
        return

    from emtrace.wrappers.minimap2 import Minimap2
    mm2 = Minimap2(num_threads=cfg.runtime.num_threads)
    mm2.build_index(db_path, mmi_path, k=mm_cfg.k, w=mm_cfg.w)
    logger.info("Host filter index ready: %s", mmi_path)


def prebuild_host_filters(
    cfg,
    keys: Optional[List[str]] = None,
    *,
    force: bool = False,
    progress_cb: ProgressCallback = None,
) -> None:
    """Pre-build minimap2 indices for several host filters.

    *keys* is a list of host database keys; ``None`` or ``["hostdbs"]`` selects
    every configured host database. Per-key failures are logged and skipped so
    one bad entry doesn't abort the rest (matching the ``hostdbs`` branch of
    :func:`install_by_key`).
    """
    if not keys or [k.lower() for k in keys] == ["hostdbs"]:
        keys = cfg.list_host_databases()

    for key in keys:
        try:
            prebuild_host_index(cfg, key, force=force, progress_cb=progress_cb)
        except DatabaseAlreadyInstalledError as exc:
            logger.warning("Skipping '%s': %s", key, exc)
        except DatabaseNotConfiguredError as exc:
            logger.warning("Skipping '%s': %s", key, exc)
        except Exception as exc:
            logger.error("Failed to build index for '%s': %s", key, exc)


def install_host_database(
    key: str,
    url: str,
    db_path: Path,
    progress_cb: ProgressCallback = None,
    force: bool = False,
) -> None:
    """Compatibility shim — delegates to :func:`install_host_minimap2`."""
    install_host_minimap2(key, [url] if url and url not in ("NA", "0") else [],
                          db_path, progress_cb, force)


def install_ncbi_taxonomy(
    sqlite_db: Path,
    url: str = (
        "ftp://ftp.ncbi.nlm.nih.gov/"
        "pub/taxonomy/taxdump.tar.gz"
    ),
    progress_cb: ProgressCallback = None,
    force: bool = False,
    version: Optional[str] = None,
) -> None:
    """Download NCBI taxdump and build the SQLite taxonomy database.

    Steps:
      1. Download taxdump.tar.gz → *sqlite_db*.parent/taxdump.tar.gz.tmp
      2. Extract → nodes.dmp, names.dmp, division.dmp, merged.dmp
      3. Call ``NCBITaxonomy.build_database()``

    Parameters
    ----------
    sqlite_db : Path
        Full path to the SQLite DB file (e.g. ``/data/taxonomy/ncbi/taxonomy.db``).
        Dump files are downloaded into ``sqlite_db.parent``.
    url : str
        NCBI taxdump URL.
    progress_cb : callable, optional
        Progress callback.
    force : bool
        Overwrite existing SQLite DB if True.
    version : str, optional
        Version string for db_info table (default: today's date).
    """
    db_dir = sqlite_db.parent
    if not force and sqlite_db.exists():
        raise DatabaseAlreadyInstalledError(
            "taxonomy", sqlite_db,
        )
    db_dir.mkdir(parents=True, exist_ok=True)

    # Phase 1: download + extract dump files
    tmp = download_to_temp(url, db_dir, progress_cb)
    extract_tarball(tmp, db_dir)

    # Phase 2: build SQLite (indeterminate progress)
    if progress_cb is not None:
        progress_cb(0, -1)

    ver = version or datetime.datetime.now().strftime("v%Y%m%d")
    logger.info("Building NCBI taxonomy DB (version=%s)...", ver)

    from emtrace.models.taxonomy import NCBITaxonomy
    tax = NCBITaxonomy(sqlite_db)
    tax.build_database(
        nodes_dmp=db_dir / "nodes.dmp",
        names_dmp=db_dir / "names.dmp",
        division_dmp=db_dir / "division.dmp",
        merged_dmp=db_dir / "merged.dmp",
        virus_only=False,
        version=ver,
    )
    logger.info("NCBI taxonomy DB ready: %s", sqlite_db)


def install_silva_taxonomy(
    sqlite_db: Path,
    taxmap_path: Path,
    url: str,
    progress_cb: ProgressCallback = None,
    force: bool = False,
    tax_slv_url: str | None = None,
    tax_slv_path: Path | None = None,
) -> None:
    """Download the SILVA taxmap file and build the SQLite taxonomy database.

    Parameters
    ----------
    sqlite_db : Path
        Full path to the SQLite DB file (e.g. ``/data/taxonomy/silva/taxonomy.db``).
        The taxmap file is downloaded into ``sqlite_db.parent``.
    taxmap_path : Path
        Destination path for the decompressed taxmap TSV.
    url : str
        URL to the ``taxmap_slv_ssu_ref_*.txt.gz`` file.
    """
    db_dir = sqlite_db.parent
    if not force and sqlite_db.exists():
        raise DatabaseAlreadyInstalledError("silva_taxonomy", sqlite_db)

    if url in ("NA", "", "0"):
        raise DatabaseNotConfiguredError("silva_taxonomy")

    db_dir.mkdir(parents=True, exist_ok=True)

    # Phase 1: download + decompress taxmap
    tmp = download_to_temp(url, db_dir, progress_cb)
    if _is_plain_gz(url):
        decompress_gz(tmp, taxmap_path)
    else:
        move_file(tmp, taxmap_path)
    logger.info("SILVA taxmap installed: %s", taxmap_path)

    # Phase 1b: the rank file. SILVA declares the rank of every path prefix here.
    # From release 144 the lineage gained a `kingdom` level, so positional rank
    # assignment puts the family name in the genus column — this file is what
    # keeps the parse correct across releases.
    if tax_slv_url and tax_slv_url not in ("NA", "", "0"):
        if tax_slv_path is None:
            tax_slv_path = db_dir / "tax_slv_ssu.txt"
        tmp_rank = download_to_temp(tax_slv_url, db_dir, progress_cb)
        if _is_plain_gz(tax_slv_url):
            decompress_gz(tmp_rank, tax_slv_path)
        else:
            move_file(tmp_rank, tax_slv_path)
        logger.info("SILVA rank file installed: %s", tax_slv_path)
    else:
        tax_slv_path = None
        logger.warning(
            "No taxonomy.silva.tax_slv_url configured — ranks will be assigned "
            "positionally, which is incorrect for SILVA 144 and later."
        )

    # Phase 2: build SQLite (indeterminate progress)
    if progress_cb is not None:
        progress_cb(0, -1)

    ver = datetime.datetime.now().strftime("v%Y%m%d")
    logger.info("Building SILVA taxonomy DB (version=%s)...", ver)

    from emtrace.models.taxonomy import SILVATaxonomy
    tax = SILVATaxonomy(sqlite_db)
    tax.build_database(
        taxmap_path=taxmap_path, tax_path=tax_slv_path, version=ver,
    )
    logger.info("SILVA taxonomy DB ready: %s", sqlite_db)


# ------------------------------------------------------------------
# Dispatcher — resolves key to database entry and installs
# ------------------------------------------------------------------

def install_by_key(
    cfg,
    key: str,
    progress_cb: ProgressCallback = None,
    force: bool = False,
) -> None:
    """Install a database identified by *key*.

    Accepts:
    - Annotation database key: ``"minimap.nt.abv"``, ``"blastn.nt.vi"``
    - Host database key, latinName, or commonName: ``"Homo_sapiens"``
    - ``"taxonomy"`` — NCBI taxonomy
    - ``"silva"`` — SILVA taxonomy
    - ``"hostdbs"`` — all configured host genomes

    Parameters
    ----------
    cfg : EMtraceConfig
        Loaded configuration.
    key : str
        Database identifier.
    progress_cb : callable, optional
        Progress callback ``(bytes_done, bytes_total)``.
    force : bool
        Overwrite existing files.
    """
    key_lower = key.lower()

    # ---- Taxonomy ---
    if key_lower == "taxonomy":
        sqlite_db = Path(cfg.taxonomy.NCBI.db)
        url = cfg.taxonomy.NCBI.url or (
            "ftp://ftp.ncbi.nlm.nih.gov/"
            "pub/taxonomy/taxdump.tar.gz"
        )
        install_ncbi_taxonomy(
            sqlite_db, url=url,
            progress_cb=progress_cb, force=force,
        )
        return

    if key_lower in ("silva", "silva_taxonomy"):
        silva_cfg = cfg.taxonomy.SILVA
        if not silva_cfg.db:
            raise DatabaseNotConfiguredError(
                "silva_taxonomy (set taxonomy.silva.db in your config)"
            )
        sqlite_db = Path(silva_cfg.db)
        # Derive taxmap path from config or use a default next to the SQLite DB
        if silva_cfg.taxmap_slv:
            taxmap_path = Path(silva_cfg.taxmap_slv)
        else:
            taxmap_path = sqlite_db.parent / "taxmap_slv_ssu_ref.txt"
        url = silva_cfg.url or ""
        tax_slv_path = (
            Path(silva_cfg.tax_slv) if silva_cfg.tax_slv
            else sqlite_db.parent / "tax_slv_ssu.txt"
        )
        install_silva_taxonomy(
            sqlite_db=sqlite_db,
            taxmap_path=taxmap_path,
            url=url,
            progress_cb=progress_cb,
            force=force,
            tax_slv_url=silva_cfg.tax_slv_url,
            tax_slv_path=tax_slv_path,
        )
        return

    # ---- All host genomes ---
    if key_lower == "hostdbs":
        for hkey in cfg.list_host_databases():
            try:
                _install_host_by_key(
                    cfg, hkey, progress_cb, force,
                )
            except (
                DatabaseAlreadyInstalledError,
                DatabaseNotConfiguredError,
            ) as exc:
                logger.warning(
                    "Skipping '%s': %s", hkey, exc,
                )
            except Exception as exc:
                logger.error(
                    "Failed to install '%s': %s", hkey, exc,
                )
        return

    # ---- Annotation databases ---
    if key in cfg.databases:
        url = cfg.get_db_url(key) or ""
        db_path = cfg.resolve_db_path(key)
        if db_path is None:
            raise DatabaseInstallError(
                f"Cannot resolve path for '{key}': "
                "EMTRACE_DB_DIR environment variable is not set. "
                "Set it in the Database Manager before installing."
            )
        install_ann_database(
            key=key,
            url=url,
            db_path=db_path,
            progress_cb=progress_cb,
            force=force,
        )
        return
    # Search by name
    for k, entry in cfg.databases.items():
        if entry.name == key:
            db_path = cfg.resolve_db_path(k)
            if db_path is None:
                raise DatabaseNotConfiguredError(k)
            install_ann_database(
                key=k,
                url=entry.url or "",
                db_path=db_path,
                progress_cb=progress_cb,
                force=force,
            )
            return

    # ---- Host database (by key, latinName, commonName) ---
    matched_host = _find_host_key(cfg, key)
    if matched_host:
        _install_host_by_key(
            cfg, matched_host, progress_cb, force,
        )
        return

    raise DatabaseNotConfiguredError(key)


def _find_host_key(cfg, key: str) -> Optional[str]:
    """Return config key matching host DB key or name."""
    if key in cfg.host_databases:
        return key
    for k, entry in cfg.host_databases.items():
        if entry.name == key:
            return k
    return None


def _install_host_by_key(
    cfg,
    key: str,
    progress_cb: ProgressCallback,
    force: bool,
) -> None:
    entry = cfg.host_databases[key]
    db_path = Path(os.path.expandvars(os.path.expanduser(entry.db)))
    install_host_minimap2(
        key=key,
        urls=entry.url_list,
        db_path=db_path,
        progress_cb=progress_cb,
        force=force,
    )


# ------------------------------------------------------------------
# Status helpers for CLI --list
# ------------------------------------------------------------------

def taxonomy_status(cfg) -> List[dict]:
    """Return status rows for NCBI and SILVA taxonomy databases.

    Each row dict has keys: key, name, url, path, installed.
    """
    rows = []

    # NCBI
    ncbi = cfg.taxonomy.NCBI
    sqlite_n = Path(ncbi.db) if ncbi.db else None
    rows.append({
        "key": "taxonomy",
        "name": ncbi.name or "NCBI taxonomy",
        "url": ncbi.url or "",
        "path": str(sqlite_n) if sqlite_n else "",
        "installed": sqlite_n.exists() if sqlite_n else False,
    })

    # SILVA
    silva = cfg.taxonomy.SILVA
    sqlite_s = Path(silva.db) if silva.db else None
    rows.append({
        "key": "silva",
        "name": silva.name or "SILVA taxonomy (SSU)",
        "url": silva.url or "",
        "path": str(sqlite_s) if sqlite_s else "",
        "installed": sqlite_s.exists() if sqlite_s else False,
    })

    return rows


def database_status(cfg) -> Tuple[List[dict], List[dict]]:
    """Return (ann_rows, host_rows) with installed status.

    Each row is a dict with keys:
      key, name, search/latinName, url, installed, path, files

    Ann rows include both new-style ``databases.sources`` entries and
    legacy ``ann.databases`` entries (deduplicated by key).
    """
    ann_rows = []
    for key in sorted(cfg.databases):
        entry = cfg.databases[key]
        db_path = cfg.resolve_db_path(key)
        installed = cfg.is_db_installed(key)
        ann_rows.append({
            "key": key,
            "name": entry.name,
            "search": entry.search,
            "url": entry.url or "",
            "path": str(db_path) if db_path else "",
            "installed": installed,
        })

    host_rows = []
    for key in cfg.list_host_databases():
        entry = cfg.host_databases[key]
        db_path = cfg.resolve_db_path(key) or Path(
            os.path.expandvars(os.path.expanduser(entry.db))
        )
        installed = check_host_installed(db_path)
        urls = entry.url_list
        # url: display string for the Source URL column
        if len(urls) > 1:
            url_display = f"{len(urls)} files"
        elif urls:
            url_display = urls[0]
        else:
            url_display = ""
        host_rows.append({
            "key": key,
            "name": entry.name,
            "search": entry.search,
            "urls": urls,
            "url": url_display,
            "path": str(db_path),
            "installed": installed,
        })

    return ann_rows, host_rows
