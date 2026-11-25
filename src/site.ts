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
    "I'm an undergraduate from Taiwan studying <strong>Economics and Statistics at Williams College (Class of 2029)</strong>. I am interested in how AI can inform human decision making, including causal inference, and how will it impact growth and inequality.", 

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

// The Chinese (繁體中文) version of the site, shown at /zh/. Anything left out
// here uses the English value above. Links in `nav` are turned into their /zh/
// versions automatically.
export const siteZh: Partial<typeof site> = {
  description: 'Haley (Yu-Heng) Lin 的個人網站：總體經濟學相關的研究與寫作。',
  role: '經濟學與統計學雙主修大學生',
  affiliation: '威廉斯學院（Williams College）',
  location: '美國麻州威廉斯敦',
  quote: '所有的判斷，就其本質而言，都是統計。——C. R. Rao',
  bio: [
    '我是來自台灣的大學生，目前就讀<strong>威廉斯學院（Williams College）經濟學與統計學系（2029 級）</strong>。我關心 AI 如何協助人類做決策（包括因果推論），以及 AI 將如何影響經濟成長與不平等。',

    '我計畫攻讀總體經濟學或 AI 的博士學位。這個網站收錄了我的研究，還有一個部落格，讓我把正在讀、正在想的其他事情寫下來。',

    '閒暇時，我喜歡寫短篇小說和玩桌遊。我在《現代藝術》（Modern Art）這款桌遊上，已經連勝三年了。',
  ],
  interests: ['總體經濟學', '因果推論', 'AI 賦能科學'],
  nav: [
    { label: '研究', href: '/research/' },
    { label: '部落格', href: '/blog/' },
    { label: '履歷', href: '/cv.pdf' },
  ],
};

/** Everything in `site`, in the given language. */
export const siteFor = (lang: 'en' | 'zh'): typeof site =>
  lang === 'zh' ? { ...site, ...siteZh } : site;
