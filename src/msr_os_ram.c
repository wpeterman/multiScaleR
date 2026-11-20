/*
 * Total physical memory reported by the operating system, in bytes.
 *
 * multiScaleR compares planned automatic-differentiation work against total
 * RAM before it builds a tape. Asking the operating system directly avoids
 * starting an outside process (PowerShell on Windows), which can be slow,
 * time out on busy machines, or be blocked by policy.
 *
 * This file deliberately includes no R headers, so that windows.h cannot clash
 * with R's macros. The .Call entry point that wraps it lives in init.c.
 * Every branch returns -1 when the query fails; the R side turns that into NA.
 */

#if defined(_WIN32)

#include <windows.h>

double msr_os_total_ram(void)
{
    MEMORYSTATUSEX status;
    status.dwLength = sizeof(status);
    if (GlobalMemoryStatusEx(&status) && status.ullTotalPhys > 0)
        return (double) status.ullTotalPhys;
    return -1.0;
}

#elif defined(__APPLE__)

#include <stdint.h>
#include <sys/types.h>
#include <sys/sysctl.h>

double msr_os_total_ram(void)
{
    int64_t mem = 0;
    size_t len = sizeof(mem);
    if (sysctlbyname("hw.memsize", &mem, &len, NULL, 0) == 0 && mem > 0)
        return (double) mem;
    return -1.0;
}

#else

#include <unistd.h>

double msr_os_total_ram(void)
{
#if defined(_SC_PHYS_PAGES) && defined(_SC_PAGESIZE)
    long pages = sysconf(_SC_PHYS_PAGES);
    long page_size = sysconf(_SC_PAGESIZE);
    if (pages > 0 && page_size > 0)
        return (double) pages * (double) page_size;
#endif
    return -1.0;
}

#endif
