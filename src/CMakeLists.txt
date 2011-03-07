
file (GLOB FILES base/*.c)
foreach (FILE ${FILES})
  set (basefiles ${basefiles} ${FILE})
endforeach ()

file (GLOB FILES hilev/*.c)
foreach (FILE ${FILES})
  get_filename_component(RELFILE ${FILE} NAME) 
  if (${RELFILE} MATCHES ExodusII)
    if (ENABLE_ExodusII)
      set (hilevfiles ${hilevfiles} ${FILE})
    endif ()
  else ()
    set (hilevfiles ${hilevfiles} ${FILE})
  endif ()
endforeach ()

file (GLOB FILES misc/*.c)
foreach (FILE ${FILES})
  set (miscfiles ${miscfiles} ${FILE})
endforeach ()

file (GLOB FILES util/*.c)
foreach (FILE ${FILES})
  set (utilfiles ${utilfiles} ${FILE})
endforeach ()

if (ENABLE_PARALLEL)
  file (GLOB FILES par/*.c)
  foreach (FILE ${FILES})
    set (parfiles ${parfiles} ${FILE})
  endforeach ()
endif ()


if (ENABLE_PARALLEL)
  set(mstkfiles ${basefiles} ${hilevfiles} ${miscfiles} ${utilfiles} ${parfiles} PARENT_SCOPE)
else ()
  set(mstkfiles ${basefiles} ${hilevfiles} ${miscfiles} ${utilfiles} PARENT_SCOPE)
endif ()

