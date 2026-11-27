@vertex
fn vs_main(@location(0) idx: u32) -> @builtin(position) vec4<f32> {
    var vertex: array<vec2<f32>, 4>;
    vertex[0] = vec2<f32>(-1.0, -1.0);
    vertex[1] = vec2<f32>(1.0, -1.0);
    vertex[2] = vec2<f32>(-1.0, 1.0);
    vertex[3] = vec2<f32>(1.0, 1.0);
    return vec4<f32>(vertex[idx], 0.0, 1.0);
}

@fragment
fn fs_main(@builtin(position) position: vec4<f32>) -> @location(0) vec4<f32> {
    var env: Environment;
    env.resolution = vec2<f32>(info__.resolution);
    env.mouse = mouse__.mouse;
    env.time = info__.time;
    let coord = vec2<f32>(
        position.x,
        env.resolution.y - position.y,
    );
    return main_image(coord, env);
}
