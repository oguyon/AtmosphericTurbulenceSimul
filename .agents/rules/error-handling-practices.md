---
description: Standardized error handling, exit codes, and clean goto cleanup patterns.
---

# Error Handling Practices

All C code in `milkatmturb` must adhere to consistent, robust, and clean error handling.

## 1. Traceable Logging
- Avoid raw uncontextualized prints for errors.
- Always log error messages to `stderr` with clear context (module, line, function):
  ```c
  fprintf(stderr, "ERROR: [%s:%d] Failed to allocate screen buffer\n", __func__, __LINE__);
  ```

## 2. Standardized Return Codes
- Functions that return a status code should return `0` on success, or a non-zero / negative value
  on failure (`errno_t`).
- Functions returning pointers should return `NULL` on failure, setting `errno` or logging a
  traceable error if necessary.

## 3. Mandatory-Check Syscalls & Libraries
The return values of file operations, memory mappings, and system calls MUST be inspected:
- `open`, `close`, `read`, `write`, `fopen`, `fclose`, `mmap`, `munmap`, `ImageStreamIO_*` calls.
- If a failure is detected, log the error using `strerror(errno)` and propagate it up the call
  stack.

## 4. Single-Label Cleanup with `goto cleanup;`
When a function holds multiple resources (descriptors, memory, ImageStreamIO handles, FFTW plans)
and may fail mid-acquisition, use the single-label `goto cleanup;` pattern with reverse-order
release at the bottom of the function. Initialise each resource handle to its sentinel (`-1`,
  `NULL`)
so the cleanup block can release them safely.

```c
int load_profile_data(const char *path, size_t sz)
{
    int   rv = -1;
    int   fd = -1;
    char *buf = NULL;

    fd = open(path, O_RDONLY);
    if (fd == -1)
    {
        fprintf(stderr, "ERROR: [%s:%d] open(%s) failed: %s\n",
                __func__, __LINE__, path, strerror(errno));
        goto cleanup;
    }

    buf = malloc(sz);
    if (buf == NULL)
    {
        fprintf(stderr, "ERROR: [%s:%d] malloc failed\n", __func__, __LINE__);
        goto cleanup;
    }

    if (read(fd, buf, sz) != (ssize_t)sz)
    {
        fprintf(stderr, "ERROR: [%s:%d] read failed: %s\n",
                __func__, __LINE__, strerror(errno));
        goto cleanup;
    }

    /* Success path */
    rv = 0;

cleanup:
    if (buf != NULL)
    {
        free(buf);
    }
    if (fd != -1)
    {
        close(fd);
    }
    return rv;
}
```
