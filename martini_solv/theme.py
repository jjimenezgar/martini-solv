"""MartiniSurf-inspired visual system for MartiniSolv."""

STYLE = """
<style>
:root {
  --ms-bg: #090A0F;
  --ms-panel: #15171F;
  --ms-input: #10131B;
  --ms-elevated: #1A1D27;
  --ms-line: #303441;
  --ms-line-strong: #405063;
  --ms-text: #F4F5F7;
  --ms-muted: #9DA3AE;
  --ms-teal: #42C7D5;
  --ms-blue: #42C7D5;
  --ms-blue-light: #8FEAF2;
  --ms-pink: #FF4FA3;
  --ms-pink-hover: #FF69B2;
  --ms-primary: var(--ms-blue);
  --ms-primary-hover: #69D7E2;
}
.stApp {
  background: var(--ms-bg);
  color: var(--ms-text);
}
.block-container {
  max-width: 1440px;
  padding-top: 4.25rem;
  padding-bottom: 3rem;
}
section[data-testid="stMain"] > div {
  overflow: visible;
}
[data-testid="stSidebar"] {
  background: #07080C;
  border-right: 1px solid var(--ms-line);
}
.ms-side-title {
  color: var(--ms-text);
  font-size: 1.25rem;
  font-weight: 850;
  letter-spacing: -0.02em;
  margin-bottom: .15rem;
}
.ms-side-flow {
  margin-top: 1rem;
  color: var(--ms-muted);
  font-size: .78rem;
  line-height: 1.5;
}

.ms-home-hero {
  position: relative;
  padding: 3.2rem 3rem 2.6rem;
  border: 1px solid rgba(66,199,213,.28);
  border-radius: 22px;
  background:
    radial-gradient(circle at 86% 16%, rgba(66,199,213,.16), transparent 30%),
    linear-gradient(135deg, rgba(18,24,34,.98), rgba(9,10,15,.98));
  box-shadow: 0 24px 70px rgba(0,0,0,.28), inset 0 0 60px rgba(66,199,213,.035);
  overflow: hidden;
  margin-bottom: 1.6rem;
}
.ms-home-hero::after {
  content: "";
  position: absolute;
  inset: auto -5rem -7rem auto;
  width: 17rem;
  height: 17rem;
  border: 1px solid rgba(143,234,242,.16);
  border-radius: 50%;
}
.ms-home-kicker {
  color: var(--ms-teal);
  font-size: .76rem;
  font-weight: 850;
  letter-spacing: .16em;
  margin-bottom: .75rem;
}
.ms-home-hero h1 {
  color: var(--ms-text);
  font-size: clamp(3rem, 8vw, 6.2rem);
  line-height: .95;
  letter-spacing: -.055em;
  margin: 0;
  font-weight: 900;
}
.ms-home-lead {
  max-width: 760px;
  margin: 1.15rem 0 1.3rem;
  color: #C8CDD5;
  font-size: 1.12rem;
  line-height: 1.65;
}
.ms-home-chips {
  display: flex;
  flex-wrap: wrap;
  gap: .55rem;
  margin: 1rem 0 2.2rem;
}
.ms-home-chips span {
  padding: .48rem .72rem;
  border: 1px solid rgba(66,199,213,.22);
  border-radius: 999px;
  background: rgba(66,199,213,.06);
  color: #DCEFF2;
  font-size: .82rem;
  font-weight: 700;
}
.ms-home-author {
  display: flex;
  gap: 1rem;
  align-items: end;
  justify-content: space-between;
  padding-top: 1.25rem;
  border-top: 1px solid rgba(157,163,174,.16);
}
.ms-home-author-label {
  color: var(--ms-muted);
  font-size: .76rem;
  text-transform: uppercase;
  letter-spacing: .09em;
}
.ms-home-author-name {
  color: var(--ms-text);
  font-size: 1.02rem;
  font-weight: 800;
  margin-top: .18rem;
}
.ms-home-github {
  color: var(--ms-teal) !important;
  text-decoration: none !important;
  font-weight: 800;
  white-space: nowrap;
}
.ms-home-github:hover {
  color: var(--ms-blue-light) !important;
}
.ms-home-card {
  min-height: 190px;
  padding: 1.25rem 1.2rem;
  border: 1px solid var(--ms-line);
  border-radius: 14px;
  background: linear-gradient(180deg, rgba(26,29,39,.92), rgba(18,20,28,.92));
}
.ms-home-card-number {
  color: var(--ms-teal);
  font-size: .74rem;
  font-weight: 900;
  letter-spacing: .1em;
  margin-bottom: 1rem;
}
.ms-home-card strong {
  display: block;
  color: var(--ms-text);
  font-size: 1.02rem;
  margin-bottom: .55rem;
}
.ms-home-card span {
  display: block;
  color: var(--ms-muted);
  font-size: .9rem;
  line-height: 1.55;
}

.ms-panel-title {
  color: var(--ms-text);
  font-size: 1.08rem;
  font-weight: 830;
  line-height: 1.35;
  padding-top: .15rem;
  margin: .05rem 0 1rem;
  overflow: visible;
}
.ms-empty-preview {
  min-height: 430px;
  border: 1px dashed var(--ms-line-strong);
  border-radius: 12px;
  background: var(--ms-panel);
  display: grid;
  place-content: center;
  text-align: center;
  padding: 2rem;
}
.ms-empty-preview strong {
  color: var(--ms-text);
  font-size: 1.25rem;
}
.ms-empty-preview span {
  display: block;
  color: var(--ms-muted);
  margin-top: .5rem;
}
div[data-testid="stMetric"] {
  background: var(--ms-panel);
  border: 1px solid var(--ms-line);
  border-radius: 10px;
  padding: .75rem .85rem;
}
div[data-testid="stMetric"] label,
div[data-testid="stMetric"] [data-testid="stMetricValue"] {
  color: var(--ms-text) !important;
}
div[data-testid="stTextInput"] input,
div[data-testid="stTextArea"] textarea,
div[data-testid="stNumberInput"] input,
div[data-baseweb="select"] > div {
  background: var(--ms-input) !important;
  color: var(--ms-text) !important;
  border-color: var(--ms-line) !important;
  border-radius: 8px !important;
}
div[data-testid="stTextInput"] input:focus,
div[data-testid="stTextArea"] textarea:focus,
div[data-testid="stNumberInput"] input:focus,
div[data-baseweb="select"] > div:focus-within {
  border-color: var(--ms-teal) !important;
  box-shadow: 0 0 0 .15rem rgba(66,199,213,.12) !important;
}
div[data-testid="stNumberInput"] button {
  background: var(--ms-elevated) !important;
  color: var(--ms-muted) !important;
  border-color: var(--ms-line) !important;
}
[data-testid="stFileUploader"] section,
[data-testid="stFileUploaderDropzone"] {
  background: var(--ms-panel) !important;
  border: 1px dashed var(--ms-line-strong) !important;
  border-radius: 10px !important;
}
[data-testid="stExpander"] {
  background: var(--ms-panel);
  border: 1px solid var(--ms-line);
  border-radius: 10px;
}
[data-testid="stExpander"] summary,
[data-testid="stExpander"] summary p {
  color: var(--ms-text) !important;
  font-weight: 760;
}
.stButton > button,
.stDownloadButton > button {
  border-radius: 8px !important;
  border: 1px solid var(--ms-line) !important;
  background: var(--ms-elevated) !important;
  color: var(--ms-text) !important;
  min-height: 2.45rem;
  font-weight: 760;
}
.stButton > button:hover,
.stDownloadButton > button:hover {
  border-color: var(--ms-teal) !important;
}
.stButton > button[kind="primary"] {
  background: var(--ms-primary) !important;
  border-color: var(--ms-primary) !important;
  color: #07131C !important;
  font-weight: 850 !important;
  box-shadow: 0 8px 22px rgba(66,199,213,.18);
}
.stButton > button[kind="primary"]:hover {
  background: var(--ms-primary-hover) !important;
  border-color: var(--ms-primary-hover) !important;
}
div[data-testid="stToggle"] [role="switch"][aria-checked="true"] {
  background-color: var(--ms-teal) !important;
}
div[data-testid="stDataFrame"] {
  border: 1px solid var(--ms-line);
  border-radius: 10px;
  overflow: hidden;
}
div[data-testid="stAlert"] {
  border-radius: 10px;
  border: 1px solid var(--ms-line);
}
label, p, span {
  color: inherit;
}
@media (max-width: 1050px) {
  .block-container {
    padding-top: 4.75rem;
    padding-left: 1rem;
    padding-right: 1rem;
  }
}
</style>
"""
