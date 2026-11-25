// English / Traditional Chinese support.
//
// English pages live at /, Chinese pages at /zh/. Every component works out
// the language from the page's URL with getLang(Astro.url), so the same
// component renders both versions.
//
// Words used across the site (buttons, labels, tag names...) are in `ui` and
// `terms` below. Your own text (bio, news, abstracts) is next to the English
// in src/site.ts, src/data/news.ts and src/data/research.ts.
// Chinese blog posts are in src/content/blog-zh/ (see README).

export type Lang = 'en' | 'zh';

export const langs: Lang[] = ['en', 'zh'];

/** Value for <html lang> and hreflang. */
export const htmlLang: Record<Lang, string> = { en: 'en', zh: 'zh-Hant-TW' };

/** Locale used to format dates. */
export const dateLocale: Record<Lang, string> = { en: 'en-US', zh: 'zh-TW' };

export const getLang = (url: URL): Lang =>
  url.pathname === '/zh' || url.pathname.startsWith('/zh/') ? 'zh' : 'en';

// Links to files (PDFs, the RSS feed, images) are the same in both languages.
const isFile = (href: string) => /\.[a-z0-9]+(?:[?#].*)?$/i.test(href.split('/').pop() ?? '');

/** Turns a site link like '/blog/' into the link for `lang` ('/zh/blog/'). */
export const localize = (href: string, lang: Lang): string => {
  if (lang === 'en' || !href.startsWith('/') || href.startsWith('//') || isFile(href)) return href;
  if (href === '/zh' || href.startsWith('/zh/')) return href;
  return `/zh${href}`;
};

/** The same page in the other language. */
export const switchPath = (pathname: string, to: Lang): string => {
  if (/^\/(zh\/)?404(\.html)?\/?$/.test(pathname)) return to === 'zh' ? '/zh/' : '/';
  const bare = pathname.replace(/^\/zh(?=\/|$)/, '') || '/';
  return localize(bare, to);
};

export const ui = {
  en: {
    langName: 'English',
    switchLabel: '中文',
    switchTitle: '切換至中文',
    mainNav: 'Main',
    themeToggle: 'Switch between light and dark mode',
    pageNotFound: 'Page not found',
    about: 'About',
    news: 'News',
    recentWriting: 'Recent writing',
    interests: 'Interests',
    selected: (heading: string) => `Selected ${heading.toLowerCase()}`,
    photoOf: (name: string) => `Photo of ${name}`,
    photoPlaceholder: 'Your photo here',
    email: 'Email',
    cv: 'CV',
    supervisedBy: 'Supervised by',
    minRead: (n: number) => `${n} min read`,
    draft: 'Draft',
    updated: 'Updated',
    contents: 'Contents',
    allPosts: '← All posts',
    filter: 'Filter',
    filterBy: 'Filter by tag',
    clear: 'Clear',
    filterEmpty: 'Nothing has all of the selected tags.',
    researchTitle: 'Research',
    researchDescription: 'Papers, thesis and other research projects.',
    researchIntro: "I'm still figuring out my interest. But here's what I've got so far.",
    researchEmpty: 'No papers have all of the selected tags.',
    blogTitle: 'Blog',
    blogDescription: 'Notes and essays on economics and other things.',
    blogIntro: 'Notes and essays: things I am learning, reading and thinking about.',
    blogEmpty: 'No posts have all of the selected tags.',
    noPosts: 'No posts yet.',
    taggedTitle: (tag: string) => `Posts tagged “${tag}”`,
    taggedDescription: (tag: string) => `Blog posts tagged ${tag}.`,
    translatedNote: '',
    originalLink: '',
  },
  zh: {
    langName: '中文',
    switchLabel: 'EN',
    switchTitle: 'Switch to English',
    mainNav: '主選單',
    themeToggle: '切換淺色／深色模式',
    pageNotFound: '找不到頁面',
    about: '關於我',
    news: '最新消息',
    recentWriting: '近期文章',
    interests: '研究興趣',
    selected: (heading: string) => `精選${heading}`,
    photoOf: (name: string) => `${name}的照片`,
    photoPlaceholder: '照片',
    email: '電子郵件',
    cv: '履歷',
    supervisedBy: '指導教授：',
    minRead: (n: number) => `閱讀時間約 ${n} 分鐘`,
    draft: '草稿',
    updated: '更新於',
    contents: '目錄',
    allPosts: '← 所有文章',
    filter: '篩選',
    filterBy: '依標籤篩選',
    clear: '清除',
    filterEmpty: '沒有同時符合所選標籤的項目。',
    researchTitle: '研究',
    researchDescription: '論文、學位論文與其他研究計畫。',
    researchIntro: '我還在摸索自己的研究興趣，以下是目前的成果。',
    researchEmpty: '沒有同時符合所選標籤的論文。',
    blogTitle: '部落格',
    blogDescription: '關於經濟學與其他事物的筆記和隨筆。',
    blogIntro: '筆記與隨筆：我正在學習、閱讀與思考的事。',
    blogEmpty: '沒有同時符合所選標籤的文章。',
    noPosts: '目前還沒有文章。',
    taggedTitle: (tag: string) => `標籤「${tag}」的文章`,
    taggedDescription: (tag: string) => `標籤為「${tag}」的部落格文章。`,
    translatedNote: '本文譯自英文原文。',
    originalLink: '閱讀英文原文',
  },
} satisfies Record<Lang, Record<string, unknown>>;

export const t = (lang: Lang) => ui[lang];

// Chinese names for tags, section headings and link labels. Anything not
// listed here is shown as written.
export const terms: Record<string, string> = {
  // Research and blog tags
  'Causal Inference': '因果推論',
  'Machine Learning': '機器學習',
  'Applied Micro': '應用個體經濟學',
  'Information Theory': '資訊理論',
  Macroeconomics: '總體經濟學',
  'AI for Science': 'AI 賦能科學',
  Personal: '個人',
  'Origin Story': '起源故事',
  // Research page groups
  Research: '研究',
  'Other Writing Samples': '其他寫作範例',
  Publications: '已發表論文',
  'Working papers': '工作論文',
  Thesis: '學位論文',
  'Writing samples': '寫作範例',
  'Other writing': '其他文章',
  // Link buttons
  Blog: '部落格',
  'Writing Samples': '寫作範例',
  PDF: 'PDF',
  Poster: '海報',
  Code: '程式碼',
  Slides: '投影片',
  Data: '資料',
  CV: '履歷',
};

/** A tag, heading or label in `lang`. */
export const term = (text: string, lang: Lang) => (lang === 'zh' ? (terms[text] ?? text) : text);

// 'Spring 2026' -> '2026 年春季', 'Oct 2026' -> '2026 年 10 月'. Plain years stay as they are.
const seasonsZh: Record<string, string> = { winter: '冬季', spring: '春季', summer: '夏季', fall: '秋季', autumn: '秋季' };
const monthsEn = ['jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec'];
export const localDate = (date: string | number, lang: Lang): string => {
  const text = String(date);
  if (lang === 'en') return text;
  const year = text.match(/\d{4}/)?.[0];
  if (!year || text.trim() === year) return text;
  const lower = text.toLowerCase();
  const season = Object.keys(seasonsZh).find((name) => lower.includes(name));
  if (season) return `${year} 年${seasonsZh[season]}`;
  const month = monthsEn.findIndex((name) => lower.includes(name));
  if (month >= 0) return `${year} 年 ${month + 1} 月`;
  return text;
};
