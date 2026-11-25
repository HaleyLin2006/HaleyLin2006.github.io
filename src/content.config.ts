import { defineCollection, z } from 'astro:content';
import { glob } from 'astro/loaders';

const schema = z.object({
  title: z.string(),
  date: z.coerce.date(),
  description: z.string().optional(),
  updated: z.coerce.date().optional(),
  tags: z.array(z.string()).default([]),
  draft: z.boolean().default(false),
});

// Blog posts are Markdown files in src/content/blog/.
// The fields above are what goes between the --- lines at the top of each post.
const blog = defineCollection({
  loader: glob({ pattern: '**/*.md', base: './src/content/blog' }),
  schema,
});

// Chinese translations of posts, shown on the /zh/ site. Give a translation the
// same file name as the post in src/content/blog/ and keep its tags in English
// (they are translated for display). Posts without a translation are shown on
// the Chinese site as they are.
const blogZh = defineCollection({
  loader: glob({ pattern: '**/*.md', base: './src/content/blog-zh' }),
  schema,
});

export const collections = { blog, blogZh };
