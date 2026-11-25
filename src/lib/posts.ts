import { getCollection, type CollectionEntry } from 'astro:content';

export type Post = CollectionEntry<'blog'>;

/** All posts, newest first. Drafts show up in `npm run dev` but not on the live site. */
export async function getPosts(): Promise<Post[]> {
  const posts = await getCollection('blog', ({ data }) =>
    import.meta.env.PROD ? !data.draft : true,
  );
  return posts.sort((a, b) => b.data.date.valueOf() - a.data.date.valueOf());
}

export function readingTime(body: string | undefined): number {
  const words = (body ?? '').trim().split(/\s+/).filter(Boolean).length;
  return Math.max(1, Math.round(words / 220));
}

export function formatDate(date: Date, withYear = true): string {
  return date.toLocaleDateString('en-US', {
    year: withYear ? 'numeric' : undefined,
    month: 'short',
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
