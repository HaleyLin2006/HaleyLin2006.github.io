---
title: The consumption Euler equation in five steps
date: 2026-11-10
description: A sample post showing how math, code and footnotes look on this site.
tags: [macro, notes]
---

This is a sample post. It is here to show how LaTeX, code blocks, tables and footnotes render, so keep it as a reference or delete it.

## The problem

A household chooses consumption $c_t$ to maximise expected lifetime utility

$$
\max_{\{c_t\}} \; \mathbb{E}_0 \sum_{t=0}^{\infty} \beta^t u(c_t)
$$

subject to the budget constraint

$$
a_{t+1} = (1 + r_{t+1})\,(a_t + y_t - c_t),
$$

where $a_t$ is wealth, $y_t$ is income and $\beta \in (0,1)$ is the discount factor.

## The first-order condition

Giving up one unit of consumption today yields $1 + r_{t+1}$ units tomorrow. At the optimum the household is indifferent to that trade:

$$
u'(c_t) = \beta \, \mathbb{E}_t \big[ (1 + r_{t+1}) \, u'(c_{t+1}) \big].
$$

This is the Euler equation.[^hall]

## With CRRA utility

Let $u(c) = \dfrac{c^{1-\sigma}}{1-\sigma}$, so that $u'(c) = c^{-\sigma}$. Then

$$
1 = \beta \, \mathbb{E}_t \left[ (1 + r_{t+1}) \left( \frac{c_{t+1}}{c_t} \right)^{-\sigma} \right].
$$

Without uncertainty, taking logs and using $\ln(1+r) \approx r$ and $\rho \equiv -\ln \beta$ gives

$$
\Delta \ln c_{t+1} \approx \frac{1}{\sigma} \, (r_{t+1} - \rho).
$$

Consumption grows when the interest rate exceeds the rate of time preference, and $1/\sigma$, the elasticity of intertemporal substitution, governs how strongly.

## A quick numerical check

```python
import numpy as np

beta, sigma, r = 0.96, 2.0, 0.05
rho = -np.log(beta)

exact = (beta * (1 + r)) ** (1 / sigma) - 1   # exact growth rate
approx = (r - rho) / sigma                    # log-linear approximation

print(f"exact: {exact:.4f}, approx: {approx:.4f}")
# exact: 0.0040, approx: 0.0046
```

## Summary

| Object | Symbol | Value in the example |
| --- | --- | --- |
| Discount factor | $\beta$ | 0.96 |
| Risk aversion | $\sigma$ | 2 |
| Interest rate | $r$ | 0.05 |

> Blockquotes look like this. Use them sparingly.

[^hall]: Hall (1978) used this condition to argue that, under some assumptions, consumption should follow a random walk.
