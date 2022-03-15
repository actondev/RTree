(defmacro rtree/defn (name &rest body)
  (let ((root (projectile-project-root) ))
    `(defun ,(intern (format "rtree/%s" name))
	 ()
       (interactive)
       (let ((default-directory ,root))
	 ,@body))))

(rtree/defn compile
	    (compile "meson compile -C build # --verbose"))

(rtree/defn test
	    (let ((debug (aod.transient/flag 'rtree "--debug")))
	      (if debug
		  (s7ty/debug-weirdNox-gdb "build/test/run-tests" "")
		(test-simple-run "meson test -C build -v"))))

(defun rtree/debug-weirdNox-gdb (command args)
  (gdb-executable command) ; nox fork ;;-i=mi
  (gdb--command (format "-exec-arguments %s" args))
  (gdb--command "-exec-run --start"))

(transient-define-prefix rtree ()
  ["rtree"
   ("q" "quit" transient-quit-all)
   ("Q" "quit & restore windows" aod-do/restore :transient transient--do-quit-all)]
  ["Args"
   ("-d" "debug" "--debug")]
  ["Commands"
   ("c" "compile" rtree/compile ; :transient transient--do-call
    )
   ("t" "test" rtree/test)
   ;;
   ])

(setq-local aod-do/action #'rtree)

(comment "custom actions"
	 
	 (aod-do/register-action #'rtree)
	 )
