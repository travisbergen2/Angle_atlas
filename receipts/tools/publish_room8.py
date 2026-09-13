#!/usr/bin/env python3
"""Create/refresh the Product_edge repo (Room VIII), push the site files, enable Pages, verify live. Idempotent.
Token from GITHUB_TOKEN (GitHub PAT Ops skill). Usage: python3 publish_room8.py [--skip-pages] [--homepage home_new2.html]"""
import os, sys, json, base64, hashlib, time, urllib.request, urllib.error
API = "https://api.github.com"; OWNER = "travisbergen2"; REPO = "Product_edge"
TOK = os.environ.get("GITHUB_TOKEN") or sys.exit("GITHUB_TOKEN not set")
HERE = os.path.dirname(os.path.abspath(__file__))
def req(method, path, body=None, ok=()):
    url = path if path.startswith("http") else API + path
    data = json.dumps(body).encode() if body is not None else None
    r = urllib.request.Request(url, data=data, method=method)
    r.add_header("Authorization", "Bearer " + TOK); r.add_header("Accept", "application/vnd.github+json"); r.add_header("X-GitHub-Api-Version", "2022-11-28")
    if data is not None: r.add_header("Content-Type", "application/json")
    try:
        resp = urllib.request.urlopen(r, timeout=60); c = resp.read(); return resp.status, (json.loads(c) if c else {})
    except urllib.error.HTTPError as e:
        if e.code in ok: return e.code, None
        sys.exit(f"HTTP {e.code} on {method} {path}: {e.read().decode(errors='replace')[:400]}")
def put_file(repo, path, local, message, branch="main"):
    raw = open(local, "rb").read()
    st, cur = req("GET", f"/repos/{repo}/contents/{path}?ref={branch}", ok=(404,))
    sha = cur.get("sha") if isinstance(cur, dict) else None
    if sha and hashlib.sha1(b"blob %d\0" % len(raw) + raw).hexdigest() == sha: print(f"  = {path} unchanged"); return
    body = {"message": message, "content": base64.b64encode(raw).decode(), "branch": branch}
    if sha: body["sha"] = sha
    st, res = req("PUT", f"/repos/{repo}/contents/{path}", body); print(f"  + {path} -> commit {res['commit']['sha'][:10]} ({len(raw)} bytes)")
def main():
    skip_pages = "--skip-pages" in sys.argv
    home = sys.argv[sys.argv.index("--homepage") + 1] if "--homepage" in sys.argv else None
    st, me = req("GET", "/user"); print("auth:", me["login"]); full = f"{OWNER}/{REPO}"
    st, r = req("GET", f"/repos/{full}", ok=(404,))
    if st == 404:
        st, r = req("POST", "/user/repos", {"name": REPO, "description": "The Product's Edge — where Euler's product over the primes stops, and why that edge is the Riemann Hypothesis. Instruments, not proofs.",
                                             "homepage": f"https://fractalyouniverse.org/{REPO}/", "has_wiki": False, "has_projects": False, "auto_init": False}); print("created repo", r["full_name"]); time.sleep(2)
    else: print("repo exists:", r["full_name"])
    msg = "The Product's Edge (Room VIII, v0.1, 2026-09-10): interactive index.html, README, receipts"
    for path in ["index.html", "README.md", "receipts/product_edge_receipts.py", "receipts/product_edge_receipts.json", "receipts/euler_product_check.py"]:
        put_file(full, path, os.path.join(HERE, path), msg)
    if not skip_pages:
        st, pg = req("GET", f"/repos/{full}/pages", ok=(404,))
        if st == 404: st, pg = req("POST", f"/repos/{full}/pages", {"source": {"branch": "main", "path": "/"}}); print("pages enabled:", pg.get("html_url"))
        else: print("pages already enabled:", pg.get("html_url"), pg.get("status"))
    if home: put_file(f"{OWNER}/travisbergen2.github.io", "index.html", home, "Room VIII — The Product's Edge: gallery card, eight instruments, description")
    local_sha = hashlib.sha256(open(os.path.join(HERE, "index.html"), "rb").read()).hexdigest()
    for url in [f"https://fractalyouniverse.org/{REPO}/", f"https://{OWNER}.github.io/{REPO}/"]:
        for i in range(12):
            try:
                body = urllib.request.urlopen(urllib.request.Request(url, headers={"Cache-Control": "no-cache", "User-Agent": "curl/8"}), timeout=20).read()
                got = hashlib.sha256(body).hexdigest(); print(f"LIVE {url} {len(body)} bytes sha256 {got[:12]} match={got == local_sha}"); break
            except Exception as e:
                if i == 11: print(f"not live yet: {url}: {e}")
                time.sleep(10)
main()
