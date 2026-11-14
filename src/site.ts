// Everything about you lives in this file. Edit it first.

export const site = {
  name: 'Haley (Yu-Heng) Lin',
  // Optional second line under your name (e.g. your name in another script). '' hides it.
  nameAlt: '',
  // Shown in browser tabs, search results and the RSS feed.
  description:
    'Personal website of Haley (Yu-Heng) Lin: research and writing on macroeconomics.',

  role: 'Undergrad in Economics and Statistics',
  affiliation: 'Williams College',
  location: 'Williamstown, MA',
  email: 'hl17@williams.edu',

  // Your photo: save it in public/images/ and set e.g. '/images/me.jpg'.
  // While this is empty, the home page shows a placeholder.
  photo: '/images/me.jpg',

  // A short line shown under the links on the home page. '' hides it.
  quote: 'All judgements are, in their rationale, statistics. - C.R. Rao',
  // Each string is one paragraph on the home page.
  // Simple HTML works: <mark>purple highlight</mark>, <strong>bold</strong>,
  // <em>italic</em>, <a href="...">link</a>.
  bio: [
    "I'm an undergraduate from Taiwan studying <strong>Economics and Statistics at Williams College (Class of 2029)</strong>. I believe AI should enhance, rather than subvert, decisions making ability of the society. I have not yet settled on my direction, but I'm interested in <strong>application of AI with theoretical guarantees</strong> and <strong>its macroeconomic consequences</strong>.",

    "I plan to pursue a PhD in macroeconomics or AI. This site collects my research, along with a blog where I think out loud about whatever else I'm reading.",
    
    "In my free time, I write short stories and play board games. I'm currently on a three-year winning streak in Modern Art.",
  ],

  // Shown as a single line under the bio. Leave empty ([]) to hide.
  interests: ['Macroeconomics', 'Causal Inference', 'AI for Science'],

  // Icon links under your photo on the home page. Delete any you do not use.
  // CV and email icons are added automatically.
  // icon is one of: 'scholar', 'github', 'linkedin', 'x'
  links: [
    // { label: 'Google Scholar', href: 'https://scholar.google.com/', icon: 'scholar' },
    { label: 'GitHub', href: 'https://github.com/HaleyLin2006', icon: 'github' },
    // { label: 'LinkedIn', href: 'https://www.linkedin.com/in/your-username', icon: 'linkedin' },
    // { label: 'X', href: 'https://x.com/your-username', icon: 'x' },
  ],

  // Put your CV at public/cv.pdf and set this to '/cv.pdf' to show a download link.
  cvPdf: '/cv.pdf',

  nav: [
    { label: 'Research', href: '/research/' },
    { label: 'Blog', href: '/blog/' },
    { label: 'CV', href: '/cv.pdf' },
  ],
};
