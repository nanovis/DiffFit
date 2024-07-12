from chimerax.core.models import Surface


class ClusterSphereModel(Surface):
    def __init__(self, name, session, color, center, radius, num_triangles=1000):
        self._num_triangles = num_triangles
        Surface.__init__(self, name, session)
        from chimerax.surface import sphere_geometry2
        va, na, ta = sphere_geometry2(self._num_triangles)
        self._unit_vertices = va
        self.set_geometry(radius*va, na, ta)
        self.color = color
        from chimerax.geometry import translation
        self.position = translation(center)
        self._radius = radius
        session.models.add([self])

    def _get_radius(self):
        return self._radius
    def _set_radius(self, r):
        if r != self._radius:
            self._radius = r
            self.set_geometry(r*self._unit_vertices, self.normals, self.triangles)
    radius = property(_get_radius, _set_radius)

