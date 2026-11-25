import { getCollection, type CollectionEntry } from 'astro:content';
import { dateLocale, type Lang } from '../i18n';

export type Post = CollectionEntry<'blog'> | CollectionEntry<'blogZh'>;

const visible = ({ data }: { data: { draft: boolean } }) =>
  import.meta.env.PROD ? !data.draft : true;

/**
 * All posts, newest first. Drafts show up in `npm run dev` but not on the live site.
 * For Chinese, a post's translation in src/content/blog-zh/ replaces it.
 */
export async function getPosts(lang: Lang = 'en'): Promise<Post[]> {
  const posts: Post[] = await getCollection('blog', visible);
  if (lang === 'zh') {
    const translations = new Map((await getCollection('blogZh', visible)).map((post) => [post.id, post]));
    for (let i = 0; i < posts.length; i++) {
      posts[i] = translations.get(posts[i].id) ?? posts[i];
    }
  }
  return posts.sort((a, b) => b.data.date.valueOf() - a.data.date.valueOf());
}

/** True for a Chinese translation of an English post. */
export const isTranslation = (post: Post) => post.collection === 'blogZh';

const cjk = /[㐀-鿿豈-﫿]/g;

/** True if the text is mostly written in Chinese. */
export const isChinese = (text: string | undefined) => (text?.match(cjk)?.length ?? 0) > 0;

/** Minutes to read: about 220 English words or 400 Chinese characters a minute. */
export function readingTime(body: string | undefined): number {
  const text = body ?? '';
  const characters = text.match(cjk)?.length ?? 0;
  const words = text.replace(cjk, ' ').trim().split(/\s+/).filter(Boolean).length;
  return Math.max(1, Math.round(words / 220 + characters / 400));
}

export function formatDate(date: Date, withYear = true, lang: Lang = 'en'): string {
  return date.toLocaleDateString(dateLocale[lang], {
    year: withYear ? 'numeric' : undefined,
    month: lang === 'zh' ? 'long' : 'short',
    day: 'numeric',
    timeZone: 'UTC',
  });
}

export function tagSlug(tag: string): string {
  return tag
    .toLowerCase()
    .trim()
    .replace(/[^a-z0-9]+/g, '-')
    .replace(/^-+|-+$/g, '');
}

export function allTags(posts: Post[]): { tag: string; slug: string; count: number }[] {
  const map = new Map<string, { tag: string; slug: string; count: number }>();
  for (const post of posts) {
    for (const tag of post.data.tags) {
      const slug = tagSlug(tag);
      if (!slug) continue;
      const entry = map.get(slug) ?? { tag, slug, count: 0 };
      entry.count += 1;
      map.set(slug, entry);
    }
  }
  return [...map.values()].sort((a, b) => a.tag.localeCompare(b.tag));
}
