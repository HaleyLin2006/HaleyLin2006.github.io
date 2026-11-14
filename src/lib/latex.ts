// Turns a short piece of LaTeX-flavoured text (e.g. an abstract) into HTML.
//
// Math is rendered with KaTeX: $inline$, \(inline\), $$display$$, \[display\].
// Outside math, the common text-mode bits work too:
//   accents      \'e \`e \"o \^o \~n \H{o} \c{c} \v{s} \={a} \u{a} \.{z}
//   dashes       -- (en dash), --- (em dash)
//   quotes       ``double'' and `single'
//   spacing      ~ (non-breaking space)
//   escapes      \$ \% \& \# \_
//   styling      \emph{...} \textit{...} \textbf{...}
import katex from 'katex';

const accents: Record<string, string> = {
  "'": '́',
  '`': '̀',
  '"': '̈',
  '^': '̂',
  '~': '̃',
  '=': '̄',
  '.': '̇',
  H: '̋',
  c: '̧',
  v: '̌',
  u: '̆',
};

const escapeHtml = (text: string) =>
  text.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');

const textMode = (text: string) => {
  const plain = text
    // \H{o}, \'{e}, \'e, \c c ...
    .replace(/\\(['`"^~=.]|[Hcvu](?![a-zA-Z]))\s*(?:\{([a-zA-Z]|\\i)\}|([a-zA-Z]))/g, (_, mark, braced, bare) => {
      const letter = (braced ?? bare) === '\\i' ? 'ı' : braced ?? bare;
      return (letter + accents[mark]).normalize('NFC');
    })
    .replace(/---/g, '—')
    .replace(/--/g, '–')
    .replace(/``/g, '“')
    .replace(/''/g, '”')
    .replace(/`/g, '‘')
    .replace(/(?<!\\)~/g, ' ')
    .replace(/\\([$%&#_])/g, '$1');

  return escapeHtml(plain)
    .replace(/\\(?:emph|textit)\{([^{}]*)\}/g, '<em>$1</em>')
    .replace(/\\textbf\{([^{}]*)\}/g, '<strong>$1</strong>');
};

const math = (tex: string, displayMode: boolean) =>
  katex.renderToString(tex, { displayMode, throwOnError: false });

// $$...$$ | \[...\] | \(...\) | $...$ (a \$ is a literal dollar sign, not math)
const mathPattern =
  /\$\$([\s\S]+?)\$\$|\\\[([\s\S]+?)\\\]|\\\(([\s\S]+?)\\\)|(?<!\\)\$((?:\\.|[^$\\])+?)\$/g;

export const latexToHtml = (source: string) => {
  let html = '';
  let last = 0;
  for (const match of source.matchAll(mathPattern)) {
    const [whole, display1, display2, inline1, inline2] = match;
    html += textMode(source.slice(last, match.index));
    html += math(display1 ?? display2 ?? inline1 ?? inline2, display1 !== undefined || display2 !== undefined);
    last = match.index + whole.length;
  }
  return html + textMode(source.slice(last));
};
