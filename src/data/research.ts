// Your papers and projects. Add an entry here and it appears on the Research
// page; set `selected: true` to also show it on the home page.

export type PaperType = 'thesis' | 'working' | 'publication' | 'writing_sample' | 'other';

export interface Paper {
  title: string;
  authors: string[]; // include yourself; your name is not specially styled
  supervisors?: string[]; // optional; shown in italics as "Supervised by ..." under the authors
  date: string | number; // shown as written, e.g. 2024, 'Spring 2024', 'Oct 2026'. Sorted by its year.
  type: PaperType;
  venue?: string; // publisher, journal, conference or status, e.g. 'AISTATS 2027', 'Preprint'
  tags?: string[]; // short topic labels; shown as coloured pills and used by the filter
  // LaTeX works: $inline$ and $$display$$ math, plus \'e, \H{o}, --, \emph{...}.
  // Write it as String.raw`...` so backslashes are kept as typed.
  abstract?: string;
  links?: { label: string; href: string }[]; // PDF, Slides, Code, Data ...
  selected?: boolean;
}

// Colour for each tag (used for research and blog tags). Tags not listed here
// get a colour picked from their name.
// Available colours: red, orange, gold, green, teal, blue, indigo, purple, magenta.
export const tagColors: Record<string, TagColor> = {
  'Causal Inference': 'blue',
  'Machine Learning': 'indigo',
  'Applied Micro': 'orange',
  'Information Theory': 'green',
};

export type TagColor =
  | 'red' | 'orange' | 'gold' | 'green' | 'teal' | 'blue' | 'indigo' | 'purple' | 'magenta';

const palette: TagColor[] = ['red', 'orange', 'gold', 'green', 'teal', 'blue', 'indigo', 'purple', 'magenta'];

export const tagColor = (tag: string): TagColor => {
  if (tagColors[tag]) return tagColors[tag];
  let hash = 0;
  for (const char of tag) hash = (hash * 31 + char.charCodeAt(0)) >>> 0;
  return palette[hash % palette.length];
};

// Sort key for a date like 2024, 'Spring 2026' or 'Oct 27, 2026': its year,
// then the order of the season or month within that year.
const within = ['jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec'];
const seasons: Record<string, number> = { winter: 0, spring: 3, summer: 6, fall: 9, autumn: 9 };
export const dateKey = (date: string | number) => {
  const text = String(date).toLowerCase();
  const year = Number(text.match(/\d{4}/)?.[0] ?? 0);
  const season = Object.keys(seasons).find((name) => text.includes(name));
  const month = within.findIndex((name) => text.includes(name));
  const part = season ? seasons[season] : month >= 0 ? month : 11;
  return year * 12 + part;
};

export const byDate = (a: Paper, b: Paper) => dateKey(b.date) - dateKey(a.date);

// Top-level groups on the Research page (in order), each with its own sections.
// The home page uses the same groups for selected papers.
// Set `sectionHeadings: false` to list a group's papers without section labels.
export const categories: {
  heading: string;
  sectionHeadings?: boolean;
  sections: { type: PaperType; heading: string }[];
}[] = [
  {
    heading: 'Research',
    sectionHeadings: false,
    sections: [
      { type: 'publication', heading: 'Publications' },
      { type: 'working', heading: 'Working papers' },
      { type: 'thesis', heading: 'Thesis' },
    ],
  },
  {
    heading: 'Writing Samples',
    sectionHeadings: false,
    sections: [
      { type: 'writing_sample', heading: 'Writing samples' },
      { type: 'other', heading: 'Other writing' },
    ],
  },
];

// Anchor id for a group heading, e.g. 'Writing Samples' -> 'writing-samples'.
export const slug = (text: string) =>
  text.toLowerCase().replace(/[^a-z0-9]+/g, '-').replace(/(^-|-$)/g, '');

export const papers: Paper[] = [
  {
    title: 'When Does Non-Gaussianity Identify Unfaithful Causal Structure?',
    authors: ['Haley Lin', 'Bijan Mazaheri'],
    supervisors: ['Prof. Bijan Mazaheri'],
    date: 2026,
    type: 'working',
    venue: 'Under review at AISTATS 2027',
    tags: ['Causal Inference', 'Machine Learning'],
    abstract: String.raw`Causal discovery methods sit at two extremes. Linear-Gaussian models, like nonparametric ones, recover a causal DAG only up to its Markov equivalence class and rely on faithfulness to do so; fully non-Gaussian models (LiNGAM) recover the entire graph without faithfulness, but only by assuming non-Gaussianity almost everywhere. We characterize the transition between these settings under partial non-Gaussianity. Our starting observation is that a non-Gaussian disturbance at a vertex orients exactly the edges an atomic intervention on that vertex would, which lets us import the verification theory of intervention design. The transition is then governed by two nested thresholds: non-Gaussianity at a vertex cover of the covered edges suffices to orient the graph given faithfulness, and additionally covering the unshielded triples of the recovered skeleton settles every v-structure without the collider test that path cancellation corrupts, removing the need for faithfulness. To assess how demanding these thresholds are, we treat non-Gaussianity as randomly occurring and analyze random DAGs. For orientation, the expected number of covered edges to be covered is $O(1)$ in dense Erd\H{o}s--R\'enyi and preferential-attachment DAGs, where full identification occurs with probability bounded below uniformly in the number of variables $N$; in sparse Erd\H{o}s--R\'enyi DAGs it is $\Theta(N)$ and the probability tends to zero. On benchmark Bayesian networks, a small non-Gaussian fraction orients most edges, while exact identification of the largest networks, which have many covered edges, requires most vertices to be non-Gaussian.`,
    links: [
      { label: 'OpenReview', href: 'https://openreview.net/forum?id=TsJOWtssfZ&referrer=%5BAuthor%20Console%5D(%2Fgroup%3Fid%3Daistats.org%2FAISTATS%2F2027%2FConference%2FAuthors%23your-submissions)' },
      { label: 'Code', href: 'https://github.com/HaleyLin2006/When_Does_Non-Gaussianity_Identify_Unfaithful_Causal_Structure' },
    ],
    selected: true,
  },
  {
    title: 'Charitable Ladies: Female Intra-household Financial Independence in 1840s Northeast United States',
    authors: ['Haley Lin'],
    supervisors: ['Prof. Steven Nafziger'],
    date: 'Spring 2026',
    type: 'writing_sample',
    venue: 'ECON 455 term paper',
    tags: ['Applied Micro'],
    abstract:
      'Under coverture, married women in nineteenth-century America faced severe restrictions on property ownership, contracting, and control over earnings. While household-level data are limited in this period, I propose using charitable giving as a first step toward a more comprehensive assessment of female intra-household financial agency in the 1840s Northeast United States. I interpret female charitable donations through a household bargaining model and compile a new dataset from the Missionary Herald to test whether (1) wealthier and more educated states have higher female donation shares and (2) the female donation share increases after the Act. However, the results do not support the two hypotheses, since household characteristics cannot be adequately controlled by state-level covariates.',
    links: [
      { label: 'PDF', href: '/Charitable%20Ladies%20Female%20Financial%20Independence%201840s.pdf' },
      { label: 'Poster', href: '/Econ%20455%20Charitable%20Ladies.pdf' },
      { label: 'Code', href: '#' },
    ],
    selected: true,
  },
  {
    title: 'An Entropy View of the Central Limit Theorem',
    authors: ['Haley Lin'],
    supervisors: ['Dr. Pei-Yu, Tsai'],
    date: 2024,
    type: 'writing_sample',
    venue: 'High school term paper',
    tags: ['Information Theory'],
    abstract:
      'Entropy, originally as a concept in thermodynamics, was introduced for information theory by C. Shannon (1948). As a result, an entropy view of probability distribution was developed. In 1985, A. Barron gave a proof of the Central Limit Theorem based on inequalities of entropy and Fisher information. The information-theoretic approach offers a cardinal comparison of probability distributions, providing an alternative description of their convergence. These results are ensured by the convexity of entropy. This expository paper reviews the foundation of information descriptions of probability theory.',
    links: [
      { label: 'PDF', href: '/An%20Entropy%20View%20of%20the%20Central%20Limit%20Theorem.pdf' },
    ],
    selected: true,
  },
];
