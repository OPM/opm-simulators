;; <http://www.gnu.org/software/emacs/manual/html_node/emacs/Directory-Variables.html>
((c-mode . (
            ; indentation of blocks is 4 characters ...
            (c-basic-offset . 4)
            ; ... and we use SPACES, _not_ TABULATORS!
            (indent-tabs-mode . nil)
            ; if we encounter a tabulator anyway, this is interpreted
            ; as being 8 character wide
            (tab-width . 8)
            ; trailing whitespace is to be avoided, so let's make it
            ; annoying.
            (setq show-trailing-whitespace t)
            ; we don't indent inline opens, ...
            (c-set-offset 'inline-open 0)
            ; ... namespaces,
            (c-set-offset 'innamespace 0)
            ; and case labels further than the rest of the function
            ; body
            (c-set-offset 'case-label 0))))
