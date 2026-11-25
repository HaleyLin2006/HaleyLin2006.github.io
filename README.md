# Personal academic website

A small, quiet academic site built with [Astro](https://astro.build): home, research, blog (with LaTeX), and CV.

## Run it on your computer

You need [Node.js](https://nodejs.org) 20 or newer.

```sh
npm install      # once, to download dependencies
npm run dev      # preview at http://localhost:4321, reloads as you edit
npm run build    # build the final site into dist/
npm run preview  # preview the built site
```

## Where to edit things

| To change | Edit |
| --- | --- |
| Name, bio, email, profile links, nav | `src/site.ts` |
| Home-page photo, quote, icon links | `photo`, `quote`, `links` in `src/site.ts` (put the photo in `public/images/`) |
| Papers and projects | `src/data/research.ts` |
| News box on the home page | `src/data/news.ts` |
| CV (PDF) | replace `public/cv.pdf` |
| Blog posts | `src/content/blog/*.md` |
| Colours, fonts, spacing | top of `src/styles/global.css` |
| Site address (for RSS and sitemap) | `site` in `astro.config.mjs` |
| Intro line on the Research page | `researchIntro` in `src/i18n.ts` |
| Tab icon | `public/favicon.jpg` (linked in `src/layouts/Base.astro`) |

To offer a PDF of your CV, save it as `public/cv.pdf` and set `cvPdf: '/cv.pdf'` in `src/site.ts`.

## Write a blog post

Create a file in `src/content/blog/`, for example `my-first-post.md`. The file name becomes the address: `/blog/my-first-post/`.

```md
---
title: My first post
date: 2026-12-01
description: One sentence shown in search results and the RSS feed.
tags: [macro, notes]
draft: false
---

Your text here, in Markdown.
```

- `title` and `date` are required; the rest are optional.
- `draft: true` keeps a post visible in `npm run dev` and off the live site.
- Math: `$\beta$` inline, and `$$ ... $$` on their own lines for display equations.
- Code blocks, tables, footnotes (`[^1]`) and images all work. Put images in `public/images/` and reference them as `/images/figure.png`.
- A "Contents" list appears automatically on posts with three or more section headings.

`src/content/blog/euler-equation.md` shows all of these in use.

## Chinese version (中文)

The whole site also exists in Traditional Chinese at `/zh/` (e.g. `/zh/research/`). The 中文 / EN button in the header switches to the same page in the other language.

| To change the Chinese | Edit |
| --- | --- |
| Bio, role, location, quote, nav | `siteZh` in `src/site.ts` |
| News | the `zh` text of each item in `src/data/news.ts` (link to `/zh/...` pages) |
| Paper venue and abstract | `venueZh` and `abstractZh` in `src/data/research.ts` (titles stay in English) |
| Buttons, labels, page intros, tag names | `ui` and `terms` in `src/i18n.ts` |
| Blog posts | `src/content/blog-zh/*.md` |

To translate a post, add a file with the **same file name** to `src/content/blog-zh/`, with the same `date` and the same English `tags` (tags are translated for display through `terms`). The Chinese site shows the translation with a link to the original. A post with no translation, like one already written in Chinese, appears on both sites as written.

Anything without a Chinese version falls back to the English.

## Add a paper

Add an entry to the `papers` list in `src/data/research.ts`. Set `selected: true` to also show it on the home page. Put PDFs in `public/papers/` and link them as `/papers/file.pdf`.

## Publish on GitHub Pages (free)

1. Create a GitHub repository named `your-username.github.io` and push this folder to its `main` branch. Include the `package-lock.json` that `npm install` created.
2. Set `site` in `astro.config.mjs` to `https://your-username.github.io`.
3. In the repository, go to **Settings → Pages** and set **Source** to **GitHub Actions**.
4. Every push to `main` now rebuilds and publishes the site (see `.github/workflows/deploy.yml`).

Notes:

- Use the `your-username.github.io` repository name or a custom domain. If you publish from a differently named repository, the site lives under a sub-path (`/repo-name/`) and the links in this project, which start at `/`, will break.
- For a custom domain, add a file `public/CNAME` containing just the domain (e.g. `yourname.com`), set `site` to match, and follow GitHub's DNS instructions under Settings → Pages.

## Project layout

```
src/
  site.ts              your details
  data/                research.ts, news.ts
  content/blog/        posts (Markdown)
  content/blog-zh/     Chinese translations of posts
  i18n.ts              English / Chinese words and link helpers
  views/               the pages' content, shared by the English and /zh/ pages
  content.config.ts    the fields a post can have
  pages/               one file per page; blog/ holds the post and tag pages; zh/ the Chinese copies
  layouts/Base.astro   shared page shell (head, header, footer)
  components/          header, footer, paper list, post list
  lib/posts.ts         sorting, dates, tags, reading time
  styles/global.css    all styling
public/                files served as-is (favicon, PDFs, images)
```
