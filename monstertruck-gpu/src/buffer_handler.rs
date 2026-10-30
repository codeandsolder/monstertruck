use crate::{
    BindingResource, Buffer, BufferBinding, BufferHandler, BufferInitDescriptor, BufferUsages,
    CommandEncoder, Device, DeviceExt, GpuError, Pod, Zeroable,
};

impl BufferHandler {
    /// Creates a buffer handler from a slice.
    #[inline]
    pub fn from_slice<T: Sized + Pod + Zeroable, A: AsRef<[T]>>(
        vec: &A,
        device: &Device,
        usage: BufferUsages,
    ) -> Self {
        let buffer = device.create_buffer_init(&BufferInitDescriptor {
            contents: bytemuck::cast_slice(vec.as_ref()),
            usage,
            label: None,
        });
        let stride = size_of::<T>() as u64;
        let size = vec.as_ref().len() as u64 * stride;
        Self {
            buffer,
            size,
            stride,
        }
    }
    /// Returns the reference of the buffer.
    #[inline]
    #[must_use]
    pub const fn buffer(&self) -> &Buffer { &self.buffer }

    /// Returns the size of the buffer.
    #[inline]
    #[must_use]
    pub const fn size(&self) -> u64 { self.size }

    /// Creates a binding resource from buffer slice.
    #[inline]
    #[must_use]
    pub const fn binding_resource(&self) -> BindingResource<'_> {
        BindingResource::Buffer(BufferBinding {
            buffer: &self.buffer,
            offset: 0,
            size: None,
        })
    }

    /// Copy the values of buffer to `dest`.
    /// # Errors
    /// Returns [`GpuError::BufferTooSmall`] when `dest` cannot hold the source buffer.
    #[inline]
    pub fn copy_buffer(&self, encoder: &mut CommandEncoder, dest: &Self) -> Result<(), GpuError> {
        if self.size > dest.size {
            return Err(GpuError::BufferTooSmall {
                source_size: self.size,
                destination_size: dest.size,
            });
        }
        encoder.copy_buffer_to_buffer(&self.buffer, 0, &dest.buffer, 0, self.size);
        Ok(())
    }
}
