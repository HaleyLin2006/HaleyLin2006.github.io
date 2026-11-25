// Short updates shown in the News box on the home page, newest first
// (they are sorted by date automatically, so the order here does not matter).
//
// date: 'YYYY-MM', e.g. '2026-10'. Shown as "Oct 2026".
// text: one sentence. Simple HTML is allowed: <strong>, <em>, <a href="...">.
// links: optional small buttons after the text, e.g. to a paper or poster.
// zh: the same sentence in Chinese for the /zh/ site (link to /zh/... pages there).
//     Without it, the Chinese site shows the English text.

export interface NewsItem {
  date: string;
  text: string;
  zh?: string;
  links?: { label: string; href: string }[];
}

// Replace these examples with your own news.
export const news: NewsItem[] = [
  {
    date: '2026-11',
    text: 'New Post! <a href="/blog/haley-twelve-favorite-questions/">Haley’s Twelve Favorite Questions</a>',
    zh: '新文章！<a href="/zh/blog/haley-twelve-favorite-questions/">Haley 最喜歡的十二個問題</a>',
    links: [{ label: 'Blog', href: '/blog/' }],
  },
  {
    date: '2026-11',
    text: 'New Post! <a href="/blog/mun/">MUN Background Guide Archive</a>',
    zh: '新文章！<a href="/zh/blog/mun/">模擬聯合國背景文件彙整</a>',
    links: [{ label: 'Blog', href: '/blog/' }],
  },
  {
    date: '2026-11',
    text: 'The website is now live! I will be posting my research and writing here.',
    zh: '網站正式上線！我會在這裡分享我的研究與文章。',
  },
  {
    date: '2026-10',
    text: 'Submitted my first ever paper to <em>AISTATS 2027</em>!',
    zh: '投出了人生第一篇論文到 <em>AISTATS 2027</em>！',
    links: [{ label: 'Research', href: '/research/#research' }],
  },
  {
    date: '2026-05',
    text: 'Finished my ECON 455 term paper on women’s charitable giving in the 1840s.',
    zh: '完成了 ECON 455 的期末論文，主題是 1840 年代女性的慈善捐款。',
    links: [{ label: 'Writing Samples', href: '/research/#writing-samples' }],
  },
];

const months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'];

export const newsDate = (date: string, lang: 'en' | 'zh' = 'en') => {
  const [year, month] = date.split('-');
  if (lang === 'zh') return { month: `${Number(month)}月`, year: `${year}年` };
  return { month: months[Number(month) - 1] ?? '', year };
};

export const sortedNews = () => [...news].sort((a, b) => b.date.localeCompare(a.date));
