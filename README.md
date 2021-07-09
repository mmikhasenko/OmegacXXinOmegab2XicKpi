# An example of HEP analysis (modeling particle spectra)

The purpose of the repository is to deponstrate application of Julia to typical problems in High-Energy physics.
Particulatly, building parametric models, and fitting them to the data using the unbinned likelihood method.

## Physics publication "Observation of excited Ω0c baryons in Ω−b→Ξ+cK−π− decays"
 - LHCb-PAPER-2021-012
 - Arxiv: [2107.03419](https://arxiv.org/abs/2107.03419)
 - Sumbitted to PRD Letters

## Analysis pipeline
The analysis pipeline includes:
  1. Preselection of decay candidates
  2. Crude selection (strippinig)
  3. Fine selection using multivarite analysis (ROOT::TMVA)
  4. (*) Efficiency parametrization
  5. (*) Fit of the particle spectra
     - (*) m(Ξc+K−π−) spectrum
     - (*) m(Ξc+K−) spectrum
  6. A(*) ngular analysis of Ωc−** → Ξc+ K− decay

The steps marked with the (*) are presented in the repository,
the rest of the code is preserved internally by the LHCb collaboration

## Workflow

### Pipelines
The code contains a set of scripts (see [scripts/pipeline.jl](scripts/pipeline.jl))
and a settings file that controles the constant variable.
The code can be run in a single `Julia` session by executing the steps of pipeline.jl

[DVC utils](https://dvc.org/) are used to run the whole analysis (`dvc repro`),
the pipeline and dependences are described in [dvc.yaml](./dvc.yaml).

### Systematic studies
The systematic studies test of the assumptions made in the analysis by modifying the parameters/settings/code if needed and rerunning the whole analysis.
The changes are implemented as git patches (`git diff > patches/systematics1.patch`).
Once the default analysis is complete (see [results/default](results/default) folder and [plots/default](plots/default) folder),
a script [scripts/systematics.jl](scripts/systematics.jl) can be run.
It rerun the analysis for every patch and stores the result in a separate folder
```julia
for (name, patch) in list_of_systematics
    gitreset()
    apply_changes(patch)
    rerun_the_analysis()
    store_results_separately(name)
end
```

## Plotting
Using `Plots.jl` with [`GR`](https://github.com/jheinen/GR.jl) and [`PGFPlotX`](https://github.com/KristofferC/PGFPlotsX.jl) backends.
The figures are created as a `@recipe`. Many tricks are implemented to match LHCb standards. It includes:
 - Histograms with Poisson error bars
 - Frame style is box with ticks also on the the left and top sides of the plot
 - Alignment of the position of the x(y) axis legend on the most right (top) side
 - Annotation with relative coordinates
 - Order in the legend that differ to the order of plotting (data)
 

