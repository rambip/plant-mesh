import marimo

__generated_with = "0.20.4"
app = marimo.App(width="medium")


@app.cell
def _():
    from tubulin import DemoMesh

    return (DemoMesh,)


@app.cell
def _(DemoMesh):
    mesh = DemoMesh()
    mesh
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
