"""L2: resolve every DOI in the repo against the DOI content-negotiation API.

Run: python -I check_dois.py <dois.txt> <out.tsv>
Writes TSV: doi, status, title (truncated).
"""
import json, sys, time, urllib.error, urllib.request

dois = [l.strip() for l in open(sys.argv[1], encoding="utf-8") if l.strip()]
rows = []
for i, doi in enumerate(dois, 1):
    url = "https://doi.org/" + urllib.parse.quote(doi, safe="/:._-()<>[]")
    req = urllib.request.Request(url, headers={
        "Accept": "application/vnd.citationstyles.csl+json",
        "User-Agent": "q-diapedesis-sim DOI audit (mailto:cnssarc@gmail.com)",
    })
    status, title = "ERROR", ""
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req, timeout=30) as r:
                meta = json.load(r)
            status = "OK"
            title = (meta.get("title") or "")
            if isinstance(title, list):
                title = title[0] if title else ""
            title = " ".join(str(title).split())[:120]
            break
        except urllib.error.HTTPError as e:
            status = f"HTTP_{e.code}"
            if e.code in (404, 400):
                break
            time.sleep(2)
        except Exception as e:
            status = f"ERR_{type(e).__name__}"
            time.sleep(2)
    rows.append((doi, status, title))
    print(f"[{i:>3}/{len(dois)}] {status:<9} {doi}", flush=True)

with open(sys.argv[2], "w", encoding="utf-8") as f:
    f.write("doi\tstatus\ttitle\n")
    for r in rows:
        f.write("\t".join(r) + "\n")

bad = [r for r in rows if r[1] != "OK"]
print(f"\n=== {len(rows) - len(bad)}/{len(rows)} resolved; {len(bad)} did NOT ===")
for d, s, _ in bad:
    print(f"  {s:<9} {d}")
