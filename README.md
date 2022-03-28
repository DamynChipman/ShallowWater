# ShallowWater

## A simple shallow water equation simulation in MATLAB.

## Overview

The Shallow Water Equations (SWE) model depth averaged flow over large distances. The word "shallow" comes from the ratio of distance traveled over depth. The SWE are used to model physical waves like tsunamis or debris flow.

This repo contains a single main file, a `Mesh` class file, and a few function files, all written in MATLAB.

![Alt Text](animation.gif)

## The Model

[Wikipedia: Shallow Water Equations](https://en.wikipedia.org/wiki/Shallow_water_equations)

## Installation

Simple clone and run:

```bash
git clone https://github.com/damynchipman/ShallowWater.git
```

## Usage

Open in MATLAB and run the `main.m` file, or via command line:

```bash
matlab main.m
```

By default, the code generates a 2D, unstructured triangle mesh on a rectangular domain, builds the initial condition (combination of a bunch of Gaussians), and steps in time until the final time. The boundary conditions are all reflective.