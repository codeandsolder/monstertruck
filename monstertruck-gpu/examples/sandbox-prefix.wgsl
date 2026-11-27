struct SceneInfo {
    background_color: vec4<f32>,
    resolution: vec2<u32>,
    time: f32,
    nlights: u32,
}

struct Mouse {
    mouse: vec4<f32>,
}

@group(0)
@binding(2)
var<uniform> info__: SceneInfo;

@group(1)
@binding(0)
var<uniform> mouse__: Mouse;

struct Environment {
    resolution: vec2<f32>,
    mouse: vec4<f32>,
    time: f32,
}
