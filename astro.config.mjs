import { defineConfig } from 'astro/config';
import sitemap from '@astrojs/sitemap';
import remarkMath from 'remark-math';
import rehypeKatex from 'rehype-katex';

export default defineConfig({
  // Change this to your real address, e.g. 'https://your-username.github.io'
  // or your custom domain. It is used for the sitemap, RSS feed and share links.
  site: 'https://haleylin2006.github.io',

  integrations: [
    sitemap({
      // English pages at /, Chinese (Traditional) pages at /zh/.
      i18n: { defaultLocale: 'en', locales: { en: 'en', zh: 'zh-Hant-TW' } },
    }),
  ],

  markdown: {
    // LaTeX in posts: $inline$ and $$display$$
    remarkPlugins: [remarkMath],
    rehypePlugins: [rehypeKatex],
    shikiConfig: {
      themes: { light: 'github-light', dark: 'github-dark' },
    },
  },
});
