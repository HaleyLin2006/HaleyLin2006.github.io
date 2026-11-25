import rss from '@astrojs/rss';
import type { APIContext } from 'astro';
import { siteFor } from '../../site';
import { getPosts } from '../../lib/posts';

const site = siteFor('zh');

export async function GET(context: APIContext) {
  const posts = await getPosts('zh');
  return rss({
    title: `${site.name} · 部落格`,
    description: site.description,
    site: context.site!,
    customData: '<language>zh-TW</language>',
    items: posts.map((post) => ({
      title: post.data.title,
      description: post.data.description,
      pubDate: post.data.date,
      link: `/zh/blog/${post.id}/`,
    })),
  });
}
