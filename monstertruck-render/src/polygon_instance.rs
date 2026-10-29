use crate::{
    Arc, AttrVertex, BindGroup, BindGroupLayout, BindingResource, BlendState, BufferAddress,
    BufferHandler, ColorTargetState, ColorWrites, CompareFunction, DepthBiasState,
    DepthStencilState, Device, DeviceHandler, Face, FilterMode, FragmentState, FrontFace,
    MipmapFilterMode, MultisampleState, PipelineCompilationOptions, PipelineLayout,
    PolygonInstance, PolygonMode, PolygonState, PrimitiveState, PrimitiveTopology, RenderId,
    RenderPipeline, RenderPipelineDescriptor, Rendered, SamplerDescriptor, SceneDescriptor,
    StencilState, Texture, TextureFormat, TextureViewDescriptor, VertexAttribute,
    VertexBufferLayout, VertexFormat, VertexState, VertexStepMode, bind_group_util, impl_render_id,
};

impl PolygonInstance {
    /// Clone the instance as another drawn element.
    #[inline]
    #[must_use]
    pub fn clone_instance(&self) -> Self {
        Self {
            polygon: self.polygon.clone(),
            state: self.state.clone(),
            shaders: self.shaders.clone(),
            id: RenderId::generate(),
        }
    }
    /// Returns a reference to the instance descriptor.
    #[inline]
    #[must_use]
    pub const fn instance_state(&self) -> &PolygonState { &self.state }
    /// Returns the mutable reference to instance descriptor.
    #[inline]
    pub const fn instance_state_mut(&mut self) -> &mut PolygonState { &mut self.state }

    /// swap vertex buffers and index buffers
    #[inline]
    pub const fn swap_vertex(&mut self, other: &mut Self) {
        std::mem::swap(&mut self.polygon, &mut other.polygon);
    }

    #[inline]
    fn non_textured_bdl(device: &Device) -> BindGroupLayout {
        let entries = [
            PolygonState::matrix_bgl_entry(),
            PolygonState::material_bgl_entry(),
        ];
        bind_group_util::create_bind_group_layout(device, &entries)
    }

    #[inline]
    fn textured_bdl(device: &Device) -> BindGroupLayout {
        bind_group_util::create_bind_group_layout(
            device,
            &[
                PolygonState::matrix_bgl_entry(),
                PolygonState::material_bgl_entry(),
                PolygonState::textureview_bgl_entry(),
                PolygonState::sampler_bgl_entry(),
            ],
        )
    }

    #[inline]
    fn non_textured_bg(&self, device: &Device, layout: &BindGroupLayout) -> BindGroup {
        bind_group_util::create_bind_group(
            device,
            layout,
            vec![
                self.state.matrix_buffer(device).binding_resource(),
                self.state.material.buffer(device).binding_resource(),
            ],
        )
    }
    #[inline]
    fn textured_bg(
        &self,
        device: &Device,
        layout: &BindGroupLayout,
        texture: &Texture,
    ) -> BindGroup {
        let view = texture.create_view(&TextureViewDescriptor::default());
        let sampler = device.create_sampler(&SamplerDescriptor {
            mag_filter: FilterMode::Linear,
            min_filter: FilterMode::Nearest,
            mipmap_filter: MipmapFilterMode::Nearest,
            lod_min_clamp: 0.0,
            lod_max_clamp: 100.0,
            ..Default::default()
        });
        bind_group_util::create_bind_group(
            device,
            layout,
            vec![
                self.state.matrix_buffer(device).binding_resource(),
                self.state.material.buffer(device).binding_resource(),
                BindingResource::TextureView(&view),
                BindingResource::Sampler(&sampler),
            ],
        )
    }
}

impl Rendered for PolygonInstance {
    impl_render_id!(id);

    #[inline]
    fn vertex_buffer(&self, _: &DeviceHandler) -> (Arc<BufferHandler>, Option<Arc<BufferHandler>>) {
        let polygon = self.polygon.clone();
        (polygon.0, Some(polygon.1))
    }
    #[inline]
    fn bind_group_layout(&self, device_handler: &DeviceHandler) -> Arc<BindGroupLayout> {
        Arc::new(if self.state.texture.is_some() {
            Self::textured_bdl(device_handler.device())
        } else {
            Self::non_textured_bdl(device_handler.device())
        })
    }
    #[inline]
    fn bind_group(
        &self,
        device_handler: &DeviceHandler,
        layout: &BindGroupLayout,
    ) -> Arc<BindGroup> {
        Arc::new(self.state.texture.as_ref().map_or_else(
            || self.non_textured_bg(device_handler.device(), layout),
            |texture| self.textured_bg(device_handler.device(), layout, texture),
        ))
    }
    #[inline]
    fn pipeline(
        &self,
        device_handler: &DeviceHandler,
        layout: &PipelineLayout,
        scene_desc: &SceneDescriptor,
    ) -> Arc<RenderPipeline> {
        let device = device_handler.device();
        let (fragment_module, fragment_entry) = if self.state.texture.is_some() {
            (
                &self.shaders.tex_fragment_module,
                self.shaders.tex_fragment_entry,
            )
        } else {
            (&self.shaders.fragment_module, self.shaders.fragment_entry)
        };
        let cull_mode = self.state.backface_culling.then_some(Face::Back);
        let blend = Some(if self.state.material.alpha_blend {
            BlendState::ALPHA_BLENDING
        } else {
            BlendState::REPLACE
        });
        let depth_stencil = if scene_desc.backend_buffer.depth_test {
            Some(DepthStencilState {
                format: TextureFormat::Depth32Float,
                depth_write_enabled: Some(true),
                depth_compare: Some(CompareFunction::Less),
                stencil: StencilState::default(),
                bias: DepthBiasState::default(),
            })
        } else {
            None
        };
        let sample_count = scene_desc.backend_buffer.sample_count;
        let pipeline = device.create_render_pipeline(&RenderPipelineDescriptor {
            layout: Some(layout),
            vertex: VertexState {
                module: &self.shaders.vertex_module,
                entry_point: Some(self.shaders.vertex_entry),
                buffers: &[VertexBufferLayout {
                    array_stride: size_of::<AttrVertex>() as BufferAddress,
                    step_mode: VertexStepMode::Vertex,
                    attributes: &[
                        VertexAttribute {
                            format: VertexFormat::Float32x3,
                            offset: 0,
                            shader_location: 0,
                        },
                        VertexAttribute {
                            format: VertexFormat::Float32x2,
                            offset: 3 * 4,
                            shader_location: 1,
                        },
                        VertexAttribute {
                            format: VertexFormat::Float32x3,
                            offset: 2 * 4 + 3 * 4,
                            shader_location: 2,
                        },
                    ],
                }],
                compilation_options: PipelineCompilationOptions::default(),
            },
            fragment: Some(FragmentState {
                module: fragment_module,
                entry_point: Some(fragment_entry),
                targets: &[Some(ColorTargetState {
                    format: scene_desc.render_texture.format,
                    blend,
                    write_mask: ColorWrites::ALL,
                })],
                compilation_options: PipelineCompilationOptions::default(),
            }),
            primitive: PrimitiveState {
                topology: PrimitiveTopology::TriangleList,
                front_face: FrontFace::Ccw,
                cull_mode,
                polygon_mode: PolygonMode::Fill,
                ..Default::default()
            },
            depth_stencil,
            multisample: MultisampleState {
                count: sample_count,
                mask: !0,
                alpha_to_coverage_enabled: sample_count > 1,
            },
            label: None,
            multiview_mask: None,
            cache: None,
        });
        Arc::new(pipeline)
    }
}
