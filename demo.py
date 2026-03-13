import marimo

__generated_with = "0.20.4"
app = marimo.App(width="medium")


@app.cell
def _():
    from tubulin import DemoMesh, Seed
    import matplotlib.pyplot as plt

    return DemoMesh, Seed


@app.cell
def _(DemoMesh):
    DemoMesh()
    return


@app.cell
def _(Seed):
    strands = Seed().grow_plant().grow_skeleton().grow_strands()
    strands.build_mesh()
    return


if __name__ == "__main__":
    app.run()
