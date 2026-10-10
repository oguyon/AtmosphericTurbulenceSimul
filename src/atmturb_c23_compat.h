// SPDX-FileCopyrightText: 2026 Olivier Guyon et al
//
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file    atmturb_c23_compat.h
 * @brief   C23 language feature shims and C11 fallbacks for milkatmturb
 */

#ifndef ATMTURB_C23_COMPAT_H
#define ATMTURB_C23_COMPAT_H

#include <assert.h>
#include <stddef.h>
#include <stdint.h>

#if defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 202311L)
#define ATMTURB_HAVE_C23 1
#else
#define ATMTURB_HAVE_C23 0
#endif

/*
 * C23 boolean keywords: bool, true, false
 */
#if !ATMTURB_HAVE_C23
#include <stdbool.h>
#endif

/*
 * C23 nullptr pointer constant
 */
#if !ATMTURB_HAVE_C23
#ifndef nullptr
#define nullptr NULL
#endif
#endif

/*
 * C23 constexpr object specifier
 */
#if ATMTURB_HAVE_C23
#define ATMTURB_CONSTEXPR constexpr
#else
#define ATMTURB_CONSTEXPR const
#endif

/*
 * C23 standard attributes with GCC/Clang fallbacks
 */
#if ATMTURB_HAVE_C23
#if defined(__has_c_attribute) && __has_c_attribute(nodiscard)
#define ATMTURB_NODISCARD [[nodiscard]]
#else
#define ATMTURB_NODISCARD
#endif

#if defined(__has_c_attribute) && __has_c_attribute(maybe_unused)
#define ATMTURB_UNUSED [[maybe_unused]]
#else
#define ATMTURB_UNUSED
#endif

#if defined(__has_c_attribute) && __has_c_attribute(fallthrough)
#define ATMTURB_FALLTHROUGH [[fallthrough]]
#else
#define ATMTURB_FALLTHROUGH ((void)0)
#endif

#elif defined(__GNUC__) || defined(__clang__)
#define ATMTURB_NODISCARD    __attribute__((warn_unused_result))
#define ATMTURB_UNUSED       __attribute__((unused))
#if defined(__GNUC__) && (__GNUC__ >= 7)
#define ATMTURB_FALLTHROUGH  __attribute__((fallthrough))
#else
#define ATMTURB_FALLTHROUGH  ((void)0)
#endif

#else
#define ATMTURB_NODISCARD
#define ATMTURB_UNUSED
#define ATMTURB_FALLTHROUGH  ((void)0)
#endif

/*
 * Single-argument static assert
 */
#if ATMTURB_HAVE_C23
#define ATMTURB_STATIC_ASSERT(cond) static_assert(cond)
#else
#define ATMTURB_STATIC_ASSERT(cond) _Static_assert((cond), #cond)
#endif

/*
 * Bit manipulation utilities (<stdbit.h> in C23, builtins in C11)
 */
#if ATMTURB_HAVE_C23 && defined(__has_include) && __has_include(<stdbit.h>)
#include <stdbit.h>

/**
 * atmturb_is_pow2_u32 - Test if a 32-bit unsigned integer is a power of two
 * @val: Value to inspect.
 *
 * Return: true if val > 0 and val is an exact power of two, false otherwise.
 */
static inline bool atmturb_is_pow2_u32(
    uint32_t val)
{
    return stdc_has_single_bit(val);
}

/**
 * atmturb_is_pow2_u64 - Test if a 64-bit unsigned integer is a power of two
 * @val: Value to inspect.
 *
 * Return: true if val > 0 and val is an exact power of two, false otherwise.
 */
static inline bool atmturb_is_pow2_u64(
    uint64_t val)
{
    return stdc_has_single_bit(val);
}

/**
 * atmturb_clz_u32 - Count leading zero bits in a 32-bit unsigned integer
 * @val: Value to inspect.
 *
 * Return: Number of leading zeros (32 if val == 0).
 */
static inline unsigned int atmturb_clz_u32(
    uint32_t val)
{
    return stdc_leading_zeros(val);
}

/**
 * atmturb_popcount_u32 - Count number of set bits in a 32-bit unsigned integer
 * @val: Value to inspect.
 *
 * Return: Number of 1-bits in val.
 */
static inline unsigned int atmturb_popcount_u32(
    uint32_t val)
{
    return stdc_count_ones(val);
}

#else // !C23 || !<stdbit.h>

/**
 * atmturb_is_pow2_u32 - Test if a 32-bit unsigned integer is a power of two
 * @val: Value to inspect.
 *
 * Return: true if val > 0 and val is an exact power of two, false otherwise.
 */
static inline bool atmturb_is_pow2_u32(
    uint32_t val)
{
    return (val > 0U) && ((val & (val - 1U)) == 0U);
}

/**
 * atmturb_is_pow2_u64 - Test if a 64-bit unsigned integer is a power of two
 * @val: Value to inspect.
 *
 * Return: true if val > 0 and val is an exact power of two, false otherwise.
 */
static inline bool atmturb_is_pow2_u64(
    uint64_t val)
{
    return (val > 0ULL) && ((val & (val - 1ULL)) == 0ULL);
}

/**
 * atmturb_clz_u32 - Count leading zero bits in a 32-bit unsigned integer
 * @val: Value to inspect.
 *
 * Return: Number of leading zeros (32 if val == 0).
 */
static inline unsigned int atmturb_clz_u32(
    uint32_t val)
{
#if defined(__GNUC__) || defined(__clang__)
    return (val == 0U) ? 32U : (unsigned int)__builtin_clz(val);
#else
    if (val == 0U)
    {
        return 32U;
    }
    unsigned int cnt = 0U;
    while ((val & 0x80000000U) == 0U)
    {
        cnt++;
        val <<= 1U;
    }
    return cnt;
#endif
}

/**
 * atmturb_popcount_u32 - Count number of set bits in a 32-bit unsigned integer
 * @val: Value to inspect.
 *
 * Return: Number of 1-bits in val.
 */
static inline unsigned int atmturb_popcount_u32(
    uint32_t val)
{
#if defined(__GNUC__) || defined(__clang__)
    return (unsigned int)__builtin_popcount(val);
#else
    unsigned int cnt = 0U;
    while (val != 0U)
    {
        val &= (val - 1U);
        cnt++;
    }
    return cnt;
#endif
}

#endif // ATMTURB_HAVE_C23 && __has_include(<stdbit.h>)

#endif // ATMTURB_C23_COMPAT_H
