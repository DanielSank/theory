# theory

Welcome to my repository of math and physics theory.
If something won't build, doesn't make sense, or if you would like to ask any questions at all about this repository, simply [file an issue](https://github.com/DanielSank/theory/issues) :D

## Download

If you're not using Git, click the green "Clone or download" button and then choose "Download ZIP". If you are using Git, you know what to do.

## Building

The LaTeX documents in this repository should build with any standard LaTeX installation, using `pdflatex`.
Image files are stored mostly as SVG, and in order to build the documents the images must first be converted to PDF.
To do this, while in the directory continaing the document you want to build, run `/bin/make-figures`.
It may be convenient to add this repository's `/bin/` directory to your `PATH`.

## Style guide

### Files

* Directories names are `AllCapsCamelCase`.

* File names are `lower_case_with_underscores`. Acronyms are fine in file names, e.g. `SNR_calculation.tex`.

### TeX


#### Equation formatting
```
\begin{equation}
  f(x) = \int_0^x \cos(y) \, dy = \sin(x) \, .
\end{equation}
```

#### Aligned equation formatting
```
begin{align}
  f(x)
  &= \int_0^x \cos(y) \, dy \\
  &= \sin(x)
  \, .
\end{align}
```

#### Labelling and referring to equations
Label equations like this:
```
\begin{equation}
  f(x) = \sinx(x) \label{eq:source.name}
\end{equation}
```
Refer to equations like this
```
Let's refer to Eq.~(\ref{eq:source.name})

```
