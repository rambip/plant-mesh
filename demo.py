import marimo

__generated_with = "0.20.4"
app = marimo.App(width="medium")


@app.cell
def _():
    from tubulin import DemoMesh

    return (DemoMesh,)


@app.cell
def _(DemoMesh):
    DemoMesh()
    return


@app.cell
def _(tubulin):
    sk = tubulin.Seed().grow_plant().grow_skeleton().grow_strands()
    return (sk,)


@app.cell
def _(sk):
    sk
    return


@app.cell
def _(sk, skeleton_to_json):
    print(_get_viewer_html(skeleton_to_json(sk)))
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
