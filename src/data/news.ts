// Short updates shown in the News box on the home page, newest first
// (they are sorted by date automatically, so the order here does not matter).
//
// date: 'YYYY-MM', e.g. '2026-10'. Shown as "Oct 2026".
// text: one sentence. Simple HTML is allowed: <strong>, <em>, <a href="...">.
// links: optional small buttons after the text, e.g. to a paper or poster.

export interface NewsItem {
  date: string;
  text: string;
  links?: { label: string; href: string }[];
}

// Replace these examples with your own news.
export const news: NewsItem[] = [
  {
    date: '2026-11',
    text: 'New Post! <a href="/blog/haley-twelve-favorite-questions/">Haley’s Twelve Favorite Questions</a>',
    links: [{ label: 'Blog', href: '/blog/' }],
  },
  {
    date: '2026-11',
    text: 'New Post! <a href="/blog/mun/">MUN Background Guide Archive</a>',
    links: [{ label: 'Blog', href: '/blog/' }],
  },
  {
    date: '2026-11',
    text: 'The website is now live! I will be posting my research and writing here.',
  },
  {
    date: '2026-10',
    text: 'Submitted my first ever paper to <em>AISTATS 2027</em>!',
    links: [{ label: 'Research', href: '/research/#research' }],
  },
  {
    date: '2026-05',
    text: 'Finished my ECON 455 term paper on women’s charitable giving in the 1840s.',
    links: [{ label: 'Writing Samples', href: '/research/#writing-samples' }],
  },
];

const months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'];

export const newsDate = (date: string) => {
  const [year, month] = date.split('-');
  return { month: months[Number(month) - 1] ?? '', year };
};

export const sortedNews = () => [...news].sort((a, b) => b.date.localeCompare(a.date));
