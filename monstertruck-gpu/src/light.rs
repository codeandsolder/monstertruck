use crate::{
    BufferHandler, BufferUsages, Device, EuclideanSpace, Light, LightInfo, LightType, Point3,
    Vector3, narrow_f64_to_f32,
};

impl Light {
    #[inline]
    pub(super) fn light_info(&self) -> LightInfo {
        LightInfo {
            position: [
                narrow_f64_to_f32(self.position.x),
                narrow_f64_to_f32(self.position.y),
                narrow_f64_to_f32(self.position.z),
                1.0,
            ],
            color: [
                narrow_f64_to_f32(self.color.x),
                narrow_f64_to_f32(self.color.y),
                narrow_f64_to_f32(self.color.z),
                1.0,
            ],
            kind: [self.kind.into(), 0, 0, 0],
        }
    }

    /// Creates a `UNIFORM` buffer of light.
    ///
    /// This method is provided only for the advanced developer utility,
    /// and not used by [`Scene`](./struct.Scene.html).
    ///
    /// # Shader Example
    /// ```glsl
    /// layout(// binding info //) uniform Light {
    ///     vec4 position;      // the position of light, position.w == 1.0
    ///     vec4 color;         // the color of light, color.w == 1.0
    ///     uvec4 light_type;   // Point => uvec4(0, 0, 0, 0), Uniform => uvec4(1, 0, 0, 0)
    /// };
    /// ```
    #[inline]
    #[must_use]
    pub fn buffer(&self, device: &Device) -> BufferHandler {
        BufferHandler::from_slice(&[self.light_info()], device, BufferUsages::UNIFORM)
    }
}

impl Default for Light {
    #[inline]
    fn default() -> Self {
        Self {
            position: Point3::origin(),
            color: Vector3::new(1.0, 1.0, 1.0),
            kind: LightType::Point,
        }
    }
}

impl From<LightType> for usize {
    #[inline]
    fn from(light_type: LightType) -> Self {
        match light_type {
            LightType::Point => 0,
            LightType::Uniform => 1,
        }
    }
}

impl From<LightType> for u32 {
    #[inline]
    fn from(light_type: LightType) -> Self {
        match light_type {
            LightType::Point => 0,
            LightType::Uniform => 1,
        }
    }
}
