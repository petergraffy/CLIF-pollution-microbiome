#!/usr/bin/env python3
"""Build only explicitly approved source files; never copy site outputs/config."""
import hashlib
import json
from pathlib import Path
import zipfile

ROOT = Path(__file__).resolve().parents[1]
BASE = [".gitignore", "config/.gitignore", "renv/.gitignore", "README.md", "LICENSE", "requirements-buddy.txt", "renv/buddy.lock",
        "docs/buddy_testing.md", "code/README.md", "config/config_template.json",
        "config/README.md", "resources/mwas/clif_intermittent_med_categories.csv",
        "utils/config.R", "utils/clif_io.R"]
SCRIPTS = ["21_prepare_acute_mwas.R", "22_cache_mwas_exposures.py", "23_run_acute_mwas.R",
           "24_summarize_acute_mwas.R", "25_audit_mwas_power_specimens.R",
           "26_cache_acs_zcta_ses.py", "27_run_federated_mwas.R", "28_pool_federated_mwas.R",
           "29_report_federated_mwas.R", "30_audit_duration_candidates.R",
           "31_shared_exposure_checks.R", "32_calibrate_mwas.R", "33_site_preflight.R",
           "34_run_buddy_site.R", "35_buddy_smoke_test.R", "36_package_buddy_source.py", "37_site_characteristics.R", "38_audit_modifier_support.R"]
UTILS = ["mwas.R", "mwas_federated.R", "mwas_poisson.R", "mwas_power.R",
         "mwas_severity.R", "mwas_targets.R", "mwas_characteristics.R"]
TESTS = ["test_mwas.R", "test_mwas_power.R", "test_mwas_windows.R", "test_mwas_federated.R",
         "test_mwas_poisson.R", "test_mwas_characteristics.R", "test_acs_ses.py", "smoke_buddy_pipeline.R"]
FILES = sorted(BASE + ["code/"+s for s in SCRIPTS] + ["utils/"+s for s in UTILS] +
               ["tests/"+s for s in TESTS])

def build():
    payload = {}
    for name in FILES:
        path = ROOT / name
        if path.is_symlink() or not path.is_file():
            raise RuntimeError(f"Missing or unsafe source file: {name}")
        if path.resolve().parent != path.parent.resolve():
            raise RuntimeError(f"Unsafe source path: {name}")
        payload[name] = path.read_bytes()
    checksums = {name: hashlib.sha256(data).hexdigest() for name, data in payload.items()}
    manifest = {"format": "clif_mwas_buddy_source_v1", "files_sha256": checksums,
                "excludes": ["clinical data", "local config", "outputs", "caches", "Git history"]}
    out = ROOT / "dist"
    out.mkdir(exist_ok=True)
    archive = out / "clif_mwas_buddy_source.zip"
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as z:
        for name, data in payload.items():
            info = zipfile.ZipInfo(name, date_time=(2026, 10, 5, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            z.writestr(info, data)
        info = zipfile.ZipInfo("BUNDLE_MANIFEST.json", date_time=(2026, 10, 5, 0, 0, 0))
        info.compress_type = zipfile.ZIP_DEFLATED
        z.writestr(info, json.dumps(manifest, indent=2)+"\n")
    checksum = hashlib.sha256(archive.read_bytes()).hexdigest()
    (out / "clif_mwas_buddy_source.zip.sha256").write_text(f"{checksum}  {archive.name}\n")
    print(f"Source-only bundle: {archive} ({len(payload)} files, SHA256 {checksum})")

if __name__ == "__main__":
    build()
