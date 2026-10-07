---
name: imagestream-internals
description: Reference for ImageStreamIO shared memory layout and semaphore synchronization in
  milkatmturb.
---

# ImageStreamIO Internals

`ImageStreamIO` provides low-latency shared-memory data mapping and semaphore signaling
used by `milkatmturb` to stream simulated wavefront phases, amplitudes, and telemetry.

## 1. Memory Mapping (`/dev/shm/`)
Each stream is backed by a POSIX shared memory file at `/dev/shm/<name>.im.shm`. It contains:
- `IMAGE_METADATA`: Dimensions (`naxis`, `size`), data type (`datatype`), write counter (`cnt0`).
- Pixel array data (`array.raw`, `array.F`, etc.).
- Read/Write semaphores (`sem_t`).

## 2. Semaphore Protocol
- The **Writer** (simulation engine) updates wavefront phase or amplitude arrays, increments
  `cnt0`, and posts to all semaphores using `ImageStreamIO_sempost(img, -1)`.
- The **Reader** (e.g. adaptive optics controller, telemetry logger) blocks on its assigned
  semaphore index using `ImageStreamIO_semwait(img, semindex)` until notified of a new frame.

## 3. Circular Buffers
When a stream has 3 axes (`naxis = 3`), the buffer behaves as a circular ring buffer for
continuous telemetry. The current slice index is written to `md->cnt1` modulo `size[2]`.
Always read the current slice index atomically before copying frame data.
