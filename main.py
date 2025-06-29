import plotly.graph_objects as go
import numpy as np

# Load orbit data
data = np.loadtxt("build/orbit.csv", delimiter=",", skiprows=1)
x, y, z = data[:, 1], data[:, 2], data[:, 3]

fig = go.Figure(layout=go.Layout(title=go.layout.Title(text="test")))

# Earth sphere
N = 200
R = 6371.0
u = np.linspace(-np.pi, np.pi, N)
v = np.linspace(0, np.pi, N)
U, V = np.meshgrid(u, v)
X = R * np.sin(V) * np.cos(U)
Y = R * np.sin(V) * np.sin(U)
Z = R * np.cos(V)

# Plot orbit
fig.add_scatter3d(
    x=x, y=y, z=z, mode="lines", line=dict(color="red", width=4), name="Orbit"
)

# Plot sphere surface
fig.add_surface(x=X, y=Y, z=Z, colorscale="blues", opacity=0.75)

fig.update_layout(
    scene=dict(
        aspectmode="data",
        xaxis_title="X (km)",
        yaxis_title="Y (km)",
        zaxis_title="Z (km)",
    ),
    title="Satellite Orbit around Earth",
)

fig.show()
