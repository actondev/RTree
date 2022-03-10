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
	    (test-simple-run "meson test -C build -v"))

(transient-define-prefix rtree ()
  ["rtree"
   ("q" "quit" transient-quit-all)
   ("Q" "quit & restore windows" aod-do/restore :transient transient--do-quit-all)]
  ["Commands"
   ("c" "compile" rtree/compile ; :transient transient--do-call
    )
   ("t" "test" rtree/test)
   ;;
   ])

(comment "custom actions"
	 (setq-local aod-do/action #'rtree)
	 (aod-do/register-action #'rtree)
	 )
