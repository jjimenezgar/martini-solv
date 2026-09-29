"""Small visual vocabulary inspired by MartiniSurf, with a teal accent."""

STYLE = """
<style>
:root { --bg: #0C0D12; --panel: #15171F; --line: #303441;
        --accent: #38C6B4; --ink: #F4F5F7; --muted: #9DA3AE; }
.stApp { background: var(--bg); color: var(--ink); }
[data-testid="stSidebar"] { background: #090A0F; border-right: 1px solid var(--line); }
.block-container { max-width: 1180px; padding-top: 2.7rem; }
.hero { padding: 1.6rem 0 1rem; border-bottom: 1px solid var(--line); margin-bottom: 1.5rem; }
.eyebrow { color: var(--accent); font-size: .76rem; letter-spacing: .12em; font-weight: 800; }
.hero h1 { color: var(--accent); font-size: clamp(2.5rem, 5vw, 4.3rem);
           letter-spacing: -.05em; line-height: 1; margin: .4rem 0; }
.hero p, .subtle { color: var(--muted); }
.panel { background: var(--panel); border: 1px solid var(--line);
         border-radius: 14px; padding: 1.25rem; margin: 1rem 0;
         box-shadow: 0 8px 22px rgba(0,0,0,.2); }
.panel h3 { margin: 0 0 .35rem; }
.pill { display: inline-block; margin-right: .4rem; padding: .3rem .65rem;
        border: 1px solid var(--line); border-radius: 999px;
        background: #1A1D27; font-size: .8rem; color: var(--accent); }
div[data-testid="stMetric"] { background: var(--panel); border: 1px solid var(--line);
                             border-radius: 10px; padding: .7rem; }
</style>
"""

