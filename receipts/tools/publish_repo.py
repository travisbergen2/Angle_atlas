#!/usr/bin/env python3
"""Create/refresh the Angle_atlas repo, push the site files, enable GitHub Pages, verify the live page.
Token from GITHUB_TOKEN (injected by the GitHub PAT Ops skill). Idempotent: safe to re-run.
Usage: python3 publish_repo.py [--skip-pages] [--homepage home_new.html]"""
import os, sys, json, base64, hashlib, time, urllib.request, urllib.error
API = "https://api.github.com"; OWNER = "travisbergen2"; REPO = "Angle_atlas"
TOK = os.environ.get("GITHUB_TOKEN") or sys.exit("GITHUB_TOKEN not set")
HERE = os.path.dirname(os.path.abspath(__file__))

def req(method, path, body=None, ok=(), accept="application/vnd.github+json"):
    url = path if path.startswith("http") else API + path
    data = json.dumps(body).encode() if body is not None else None
    r = urllib.request.Request(url, data=data, method=method)
    r.add_header("Authorization", "Bearer " + TOK); r.add_header("Accept", accept); r.add_header("X-GitHub-Api-Version", "2022-11-28")
    if data is not None: r.add_header("Content-Type", "application/json")
    try:
        resp = urllib.request.urlopen(r, timeout=60); c = resp.read()
        return resp.status, (json.loads(c) if c else {})
    except urllib.error.HTTPError as e:
        if e.code in ok: return e.code, None
        sys.exit(f"HTTP {e.code} on {method} {path}: {e.read().decode(errors='replace')[:400]}")

def put_file(repo, path, local, message, branch="main"):
    with open(local, "rb") as f: raw = f.read()
    st, cur = req("GET", f"/repos/{repo}/contents/{path}?ref={branch}", ok=(404,))
    sha = cur.get("sha") if isinstance(cur, dict) else None
    if sha and hashlib.sha1(b"blob %d\0" % len(raw) + raw).hexdigest() == sha:
        print(f"  = {path} unchanged"); return
    body = {"message": message, "content": base64.b64encode(raw).decode(), "branch": branch}
    if sha: body["sha"] = sha
    st, res = req("PUT", f"/repos/{repo}/contents/{path}", body)
    print(f"  + {path} -> commit {res['commit']['sha'][:10]} ({len(raw)} bytes)")

def main():
    skip_pages = "--skip-pages" in sys.argv
    home = sys.argv[sys.argv.index("--homepage") + 1] if "--homepage" in sys.argv else None
    st, me = req("GET", "/user"); print("auth:", me["login"])
    full = f"{OWNER}/{REPO}"
    st, r = req("GET", f"/repos/{full}", ok=(404,))
    if st == 404:
        st, r = req("POST", "/user/repos", {"name": REPO, "description": "The Angle Atlas — every load-bearing angle in the IMM Riemann program on one dial. Instruments, not proofs.",
                                             "homepage": f"https://fractalyouniverse.org/{REPO}/", "has_wiki": False, "has_projects": False, "auto_init": False})
        print("created repo", r["full_name"]); time.sleep(2)
    else:
        print("repo exists:", r["full_name"], "default branch:", r.get("default_branch"))
    msg = "The Angle Atlas — the dial (v0.1, 2026-09-10): interactive index.html, README, receipts, figure"
    for path, local in [("index.html", "index.html"), ("README.md", "README.md"), ("receipts/angle_dial.py", "receipts/angle_dial.py"),
                        ("receipts/verify_angles.py", "receipts/verify_angles.py"), ("receipts/reciprocal_check.py", "receipts/reciprocal_check.py"),
                        ("angle_dial.png", "angle_dial.png")]:
        put_file(full, path, os.path.join(HERE, local), msg)
    if not skip_pages:
        st, pg = req("GET", f"/repos/{full}/pages", ok=(404,))
        if st == 404:
            st, pg = req("POST", f"/repos/{full}/pages", {"source": {"branch": "main", "path": "/"}}); print("pages enabled:", pg.get("html_url"))
        else:
            print("pages already enabled:", pg.get("html_url"), "status:", pg.get("status"))
    if home:
        put_file(f"{OWNER}/travisbergen2.github.io", "index.html", home, "Room VII — The Angle Atlas: gallery card, count Seven, description")
    # verify live
    local_sha = hashlib.sha256(open(os.path.join(HERE, "index.html"), "rb").read()).hexdigest()
    for url in [f"https://fractalyouniverse.org/{REPO}/", f"https://{OWNER}.github.io/{REPO}/"]:
        for i in range(12):
            try:
                body = urllib.request.urlopen(urllib.request.Request(url, headers={"Cache-Control": "no-cache", "User-Agent": "curl/8"}), timeout=20).read()
                got = hashlib.sha256(body).hexdigest()
                print(f"LIVE {url} {len(body)} bytes sha256 {got[:12]} match={got == local_sha}"); break
            except Exception as e:
                if i == 11: print(f"not live yet: {url}: {e}")
                time.sleep(10)
main()
