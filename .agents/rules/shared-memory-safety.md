---
description: Safety protocols when using ImageStreamIO shared memory streams.
---

# Shared Memory Safety

When using `ImageStreamIO` shared memory buffers for real-time wavefront and telemetry streaming:

## 1. Testing Cleanup
- Clean up test streams created under `/dev/shm/` immediately after tests finish:
  ```bash
  rm -f /dev/shm/test_*.im.shm
  ```

## 2. Stream Owner Validation
- Check ownership of `/dev/shm/*.im.shm` files to prevent overwriting active streams.
- Ensure streams are created with appropriate permissions and dimensions.

## 3. Synchronization
- Always use the semaphore protocol (`ImageStreamIO_sempost`, `ImageStreamIO_semwait`) to signal
  or block until frames are ready, avoiding busy-waiting spin loops.
- Use `ImageStreamIO_sempost(&img, -1)` to post to all active semaphores.
