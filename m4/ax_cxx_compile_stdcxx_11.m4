<?xml version="1.0" encoding="utf-8"?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en-US" lang="en-US">
<!-- git web interface version 1.7.2.5, (C) 2005-2006, Kay Sievers <kay.sievers@vrfy.org>, Christian Gierke -->
<!-- git core binaries version 1.7.2.5 -->
<head>
<meta http-equiv="content-type" content="text/html; charset=utf-8"/>
<meta name="generator" content="gitweb/1.7.2.5 git/1.7.2.5"/>
<meta name="robots" content="index, nofollow"/>
<title>Savannah Git Hosting - autoconf-archive.git/summary</title>
<link rel="stylesheet" type="text/css" href="gitweb.css"/>
<link rel="alternate" title="autoconf-archive.git - log - RSS feed" href="/gitweb/?p=autoconf-archive.git;a=rss" type="application/rss+xml" />
<link rel="alternate" title="autoconf-archive.git - log - RSS feed (no merges)" href="/gitweb/?p=autoconf-archive.git;a=rss;opt=--no-merges" type="application/rss+xml" />
<link rel="alternate" title="autoconf-archive.git - log - Atom feed" href="/gitweb/?p=autoconf-archive.git;a=atom;opt=--no-merges" type="application/atom+xml" />
<link rel="alternate" title="autoconf-archive.git - log - Atom feed (no merges)" href="/gitweb/?p=autoconf-archive.git;a=atom;opt=--no-merges" type="application/atom+xml" />
<link rel="shortcut icon" href="git-favicon.png" type="image/png" />
</head>
<body>
<div class="page_header">
<a title="git homepage" href="http://git-scm.com/"><img src="git-logo.png" width="72" height="27" alt="git" class="logo"/></a><a href="/gitweb/">git@sv</a> / <a href="/gitweb/?p=autoconf-archive.git;a=summary">autoconf-archive.git</a> / summary
</div>
<form method="get" action="/gitweb/" enctype="application/x-www-form-urlencoded">
<div class="search">
<input name="p" type="hidden" value="autoconf-archive.git" />
<input name="a" type="hidden" value="search" />
<input name="h" type="hidden" value="HEAD" />
<select name="st" >
<option selected="selected" value="commit">commit</option>
<option value="grep">grep</option>
<option value="author">author</option>
<option value="committer">committer</option>
<option value="pickaxe">pickaxe</option>
</select><sup><a href="/gitweb/?p=autoconf-archive.git;a=search_help">?</a></sup> search:
<input type="text" name="s"  />
<span title="Extended regular expression"><label><input type="checkbox" name="sr" value="1" />re</label></span></div>
</form>
<div class="page_nav">
summary | <a href="/gitweb/?p=autoconf-archive.git;a=shortlog">shortlog</a> | <a href="/gitweb/?p=autoconf-archive.git;a=log">log</a> | <a href="/gitweb/?p=autoconf-archive.git;a=commit;h=5c2ada8a24cc5bafd1c8e99ff217d7c29fa3844f">commit</a> | <a href="/gitweb/?p=autoconf-archive.git;a=commitdiff;h=5c2ada8a24cc5bafd1c8e99ff217d7c29fa3844f">commitdiff</a> | <a href="/gitweb/?p=autoconf-archive.git;a=tree">tree</a><br/>
<br/>
</div>
<div class="title">&nbsp;</div>
<table class="projects_list">
<tr id="metadata_desc"><td>description</td><td>GNU Autoconf Archive</td></tr>
<tr id="metadata_owner"><td>owner</td><td></td></tr>
<tr id="metadata_lchange"><td>last change</td><td>Thu, 17 Dec 2015 21:56:45 +0000</td></tr>
<tr class="metadata_url"><td>URL</td><td>git://git.sv.gnu.org/autoconf-archive.git</td></tr>
<tr class="metadata_url"><td></td><td>http://git.savannah.gnu.org/r/autoconf-archive.git</td></tr>
<tr class="metadata_url"><td></td><td>ssh://git.sv.gnu.org/srv/git/autoconf-archive.git</td></tr>
</table>
<div class="header">
<a class="title" href="/gitweb/?p=autoconf-archive.git;a=shortlog">shortlog</a>
</div>
<table class="shortlog">
<tr class="dark">
<td title="3 weeks ago"><i>2015-12-17</i></td>
<td class="author"><a title="Search for commits authored by Peter Simons" class="list" href="/gitweb/?p=autoconf-archive.git;a=search;s=Peter+Simons;st=author">Peter Simons</a></td><td><a class="list subject" href="/gitweb/?p=autoconf-archive.git;a=commit;h=5c2ada8a24cc5bafd1c8e99ff217d7c29fa3844f">Merge pull request #52 from swofford/fix-ax_gcc_cpuid</a> <span class="refs"> <span class="head" title="heads/master"><a href="/gitweb/?p=autoconf-archive.git;a=shortlog;h=refs/heads/master">master</a></span></span></td>
<td class="link"><a href="/gitweb/?p=autoconf-archive.git;a=commit;h=5c2ada8a24cc5bafd1c8e99ff217d7c29fa3844f">commit</a> | <a href="/gitweb/?p=autoconf-archive.git;a=commitdiff;h=5c2ada8a24cc5bafd1c8e99ff217d7c29fa3844f">commitdiff</a> | <a href="/gitweb/?p=autoconf-archive.git;a=tree;h=5c2ada8a24cc5bafd1c8e99ff217d7c29fa3844f;hb=5c2ada8a24cc5bafd1c8e99ff217d7c29fa3844f">tree</a> | <a title="in format: tar.gz" href="/gitweb/?p=autoconf-archive.git;a=snapshot;h=5c2ada8a24cc5bafd1c8e99ff217d7c29fa3844f;sf=tgz">snapshot</a></td>
</tr>
<tr class="light">
<td title="3 weeks ago"><i>2015-12-17</i></td>
<td class="author"><a title="Search for commits authored by Dave Swofford" class="list" href="/gitweb/?p=autoconf-archive.git;a=search;s=Dave+Swofford;st=author">Dave Swofford</a></td><td><a class="list subject" title="fix AX_GCC_X86_CPUID/AX_GCC_X86_CPUID_COUNT for 32-bit PIC compilations" href="/gitweb/?p=autoconf-archive.git;a=commit;h=fc56521c600a5cb9330daa52a093007dfe47f081">fix AX_GCC_X86_CPUID/AX_GCC_X86_CPUID_COUNT for 32... </a></td>
<td class="link"><a href="/gitweb/?p=autoconf-archive.git;a=commit;h=fc56521c600a5cb9330daa52a093007dfe47f081">commit</a> | <a href="/gitweb/?p=autoconf-archive.git;a=commitdiff;h=fc56521c600a5cb9330daa52a093007dfe47f081">commitdiff</a> | <a href="/gitweb/?p=autoconf-archive.git;a=tree;h=fc56521c600a5cb9330daa52a093007dfe47f081;hb=fc56521c600a5cb9330daa52a093007dfe47f081">tree</a> | <a title="in format: tar.gz" href="/gitweb/?p=autoconf-archive.git;a=snapshot;h=fc56521c600a5cb9330daa52a093007dfe47f081;sf=tgz">snapshot</a></td>
</tr>
<tr class="dark">
<td title="3 weeks ago"><i>2015-12-17</i></td>
<td class="author"><a title="Search for commits authored by Peter Simons" class="list" href="/gitweb/?p=autoconf-archive.git;a=search;s=Peter+Simons;st=author">Peter Simons</a></td><td><a class="list subject" href="/gitweb/?p=autoconf-archive.git;a=commit;h=1b9f59d95f8a7050849441233f5a38d52f24a22b">Merge pull request #51 from N0NB/master</a></td>
<td class="link"><a href="/gitweb/?p=autoconf-archive.git;a=commit;h=1b9f59d95f8a7050849441233f5a38d52f24a22b">commit</a> | <a href="/gitweb/?p=autoconf-archive.git;a=commitdiff;h=1b9f59d95f8a7050849441233f5a38d52f24a22b">commitdiff</a> | <a href="/gitweb/?p=autoconf-archive.git;a=tree;h=1b9f59d95f8a7050849441233f5a38d52f24a22b;hb=1b9f59d95f8a7050849441233f5a38d52f24a22b">tree</a> | <a title="in format: tar.gz" href="/gitweb/?p=autoconf-archive.git;a=snapshot;h=1b9f59d95f8a7050849441233f5a38d52f24a22b;sf=tgz">snapshot</a></td>
</tr>
<tr class="light">
<td title="3 weeks ago"><i>2015-12-16</i></td>
<td class="author"><a title="Search for commits authored by Nate Bargmann" class="list" href="/gitweb/?p=autoconf-archive.git;a=search;s=Nate+Bargmann;st=author">Nate Bargmann</a></td><td><a class="list subject" href="/gitweb/?p=autoconf-archive.git;a=commit;h=302201bf1bb20543c5f3336b1a913d877e7642fb">Fix second pass bug in _AX_WITH_CURSES_CHECKEXTRA</a></td>
<td class="link"><a href="/gitweb/?p=autoconf-archive.git;a=commit;h=302201bf1bb20543c5f3336b1a913d877e7642fb">commit</a> | <a href="/gitweb/?p=autoconf-archive.git;a=commitdiff;h=302201bf1bb20543c5f3336b1a913d877e7642fb">commitdiff</a> | <a href="/gitweb/?p=autoconf-archive.git;a=tree;h=302201bf1bb20543c5f3336b1a913d877e7642fb;hb=302201bf1bb20543c5f3336b1a913d877e7642fb">tree</a> | <a title="in format: tar.gz" href="/gitweb/?p=autoconf-archive.git;a=snapshot;h=302201bf1bb20543c5f3336b1a913d877e7642fb;sf=tgz">snapshot</a></td>
</tr>
<tr class="dark">
<td title="6 weeks ago"><i>2015-11-23</i></td>
<td class="author"><a title="Search for commits authored by Murray Cumming" class="list" href="/gitweb/?p=autoconf-archive.git;a=search;s=Murray+Cumming;st=author">Murray Cumming</a></td><td><a class="list subject" href="/gitweb/?p=autoconf-archive.git;a=commit;h=433b0e37b2c3cd7aa6b4551a7824a49e06875fba">AX_BOOST_PYTHON: Update for the AX_PYTHON_DEVEL change.</a></td>
<td class="link"><a href="/gitweb/?p=autoconf-archive.git;a=commit;h=433b0e37b2c3cd7aa6b4551a7824a49e06875fba">commit</a> | <a href="/gitweb/?p=autoconf-archive.git;a=commitdiff;h=433b0e37b2c3cd7aa6b4551a7824a49e06875fba">commitdiff</a> | <a href="/gitweb/?p=autoconf-archive.git;a=tree;h=433b0e37b2c3cd7aa6b4551a7824a49e06875fba;hb=433b0e37b2c3cd7aa6b4551a7824a49e06875fba">tree</a> | <a title="in format: tar.gz" href="/gitweb/?p=autoconf-archive.git;a=snapshot;h=433b0e37b2c3cd7aa6b4551a7824a49e06875fba;sf=tgz">snapshot</a></td>
</tr>
<tr class="light">
<td title="6 weeks ago"><i>2015-11-23</i></td>
<td class="author"><a title="Search for commits authored by Markus Armbruster" class="list" href="/gitweb/?p=autoconf-archive.git;a=search;s=Markus+Armbruster;st=author">Markus Armbruster</a></td><td><a class="list subject" title="AX_APPEND_COMPILE_FLAGS, AX_APPEND_LINK_FLAGS: Optional INPUT arg" href="/gitweb/?p=autoconf-archive.git;a=commit;h=e3d948b0303aabb8c016c7f3d7f1d538aa5472be">AX_APPEND_COMPILE_FLAGS, AX_APPEND_LINK_FLAGS: Optional... </a></td>
<td class="link"><a href="/gitweb/?p=autoconf-archive.git;a=commit;h=e3d948b0303aabb8c016c7f3d7f1d538aa5472be">commit</a> | <a href="/gitweb/?p=autoconf-archive.git;a=commitdiff;h=e3d948b0303aabb8c016c7f3d7f1d538aa5472be">commitdiff</a> | <a href="/gitweb/?p=autoconf-archive.git;a=tree;h=e3d948b0303aabb8c016c7f3d7f1d538aa5472be;hb=e3d948b0303aabb8c016c7f3d7f1d538aa5472be">tree</a> | <a title="in format: tar.gz" href="/gitweb/?p=autoconf-archive.git;a=snapshot;h=e3d948b0303aabb8c016c7f3d7f1d538aa5472be;sf=tgz">snapshot</a></td>
</tr>
<tr class="dark">
<td title="6 weeks ago"><i>2015-11-23</i></td>
<td class="author"><a title="Search for commits authored by Moritz Klammler" class="list" href="/gitweb/?p=autoconf-archive.git;a=search;s=Moritz+Klammler;st=author">Moritz Klammler</a></td><td><a class="list subject" href="/gitweb/?p=autoconf-archive.git;a=commit;h=3a8480edb2e6784f3d8f7e6e0b4181e70c626af0">Unify macros for checking C++11/14/... compiler support.</a></td>
<td class="link"><a href="/gitweb/?p=autoconf-archive.git;a=commit;h=3a8480edb2e6784f3d8f7e6e0b4181e70c626af0">commit</a> | <a href="/gitweb/?p=autoconf-archive.git;a=commitdiff;h=3a8480edb2e6784f3d8f7e6e0b4181e70c626af0">commitdiff</a> | <a href="/gitweb/?p=autoconf-archive.git;a=tree;h=3a8480edb2e6784f3d8f7e6e0b4181e70c626af0;hb=3a8480edb2e6784f3d8f7e6e0b4181e70c626af0">tree</a> | <a title="in format: tar.gz" href="/gitweb/?p=autoconf-archive.git;a=snapshot;h=3a8480edb2e6784f3d8f7e6e0b4181e70c626af0;sf=tgz">snapshot</a></td>
</tr>
<tr class="light">
<td title="6 weeks ago"><i>2015-11-23</i></td>
<td class="author"><a title="Search for commits authored by Peter Simons" class="list" href="/gitweb/?p=autoconf-archive.git;a=search;s=Peter+Simons;st=author">Peter Simons</a></td><td><a class="list subject" title="Merge branch 'code-coverage' of https://github.com/olaf-mandel/autoconf-archive into... " href="/gitweb/?p=autoconf-archive.git;a=commit;h=bd1d2f55014ddd9061545c048d982a3c11851d90">Merge branch 'code-coverage' of https://github.com... </a></td>
<td class="link"><a href="/gitweb/?p=autoconf-archive.git;a=commit;h=bd1d2f55014ddd9061545c048d982a3c11851d90">commit</a> | <a href="/gitweb/?p=autoconf-archive.git;a=commitdiff;h=bd1d2f55014ddd9061545c048d982a3c11851d90">commitdiff</a> | <a href="/gitweb/?p=autoconf-archive.git;a=tree;h=bd1d2f55014ddd9061545c048d982a3c11851d90;hb=bd1d2f55014ddd9061545c048d982a3c11851d90">tree</a> | <a title="in format: tar.gz" href="/gitweb/?p=autoconf-archive.git;a=snapshot;h=bd1d2f55014ddd9061545c048d982a3c11851d90;sf=tgz">snapshot</a></td>
</tr>
<tr class="dark">
<td title="6 weeks ago"><i>2015-11-23</i></td>
<td class="author"><a title="Search for commits authored by Olaf Mandel" class="list" href="/gitweb/?p=autoconf-archive.git;a=search;s=Olaf+Mandel;st=author">Olaf Mandel</a></td><td><a class="list subject" href="/gitweb/?p=autoconf-archive.git;a=commit;h=9bf65e5f34ab88e68c2133d2dab8da71989ae75a">AX_PROG_DOXYGEN: allow AC_CONFIG_FILES([Doxyfile])</a></td>
<td class="link"><a href="/gitweb/?p=autoconf-archive.git;a=commit;h=9bf65e5f34ab88e68c2133d2dab8da71989ae75a">commit</a> | <a href="/gitweb/?p=autoconf-archive.git;a=commitdiff;h=9bf65e5f34ab88e68c2133d2dab8da71989ae75a">commitdiff</a> | <a href="/gitweb/?p=autoconf-archive.git;a=tree;h=9bf65e5f34ab88e68c2133d2dab8da71989ae75a;hb=9bf65e5f34ab88e68c2133d2dab8da71989ae75a">tree</a> | <a title="in format: tar.gz" href="/gitweb/?p=autoconf-archive.git;a=snapshot;h=9bf65e5f34ab88e68c2133d2dab8da71989ae75a;sf=tgz">snapshot</a></td>
</tr>
<tr class="light">
<td title="6 weeks ago"><i>2015-11-23</i></td>
<td class="author"><a title="Search for commits authored by Olaf Mandel" class="list" href="/gitweb/?p=autoconf-archive.git;a=search;s=Olaf+Mandel;st=author">Olaf Mandel</a></td><td><a class="list subject" href="/gitweb/?p=autoconf-archive.git;a=commit;h=ebfdcdeb6508237bc73a18ad1f8614292b4aa88e">AX_PROG_DOXYGEN: support multiple Doxyfiles</a></td>
<td class="link"><a href="/gitweb/?p=autoconf-archive.git;a=commit;h=ebfdcdeb6508237bc73a18ad1f8614292b4aa88e">commit</a> | <a href="/gitweb/?p=autoconf-archive.git;a=commitdiff;h=ebfdcdeb6508237bc73a18ad1f8614292b4aa88e">commitdiff</a> | <a href="/gitweb/?p=autoconf-archive.git;a=tree;h=ebfdcdeb6508237bc73a18ad1f8614292b4aa88e;hb=ebfdcdeb6508237bc73a18ad1f8614292b4aa88e">tree</a> | <a title="in format: tar.gz" href="/gitweb/?p=autoconf-archive.git;a=snapshot;h=ebfdcdeb6508237bc73a18ad1f8614292b4aa88e;sf=tgz">snapshot</a></td>
</tr>
<tr class="dark">
<td title="6 weeks ago"><i>2015-11-23</i></td>
<td class="author"><a title="Search for commits authored by Olaf Mandel" class="list" href="/gitweb/?p=autoconf-archive.git;a=search;s=Olaf+Mandel;st=author">Olaf Mandel</a></td><td><a class="list subject" href="/gitweb/?p=autoconf-archive.git;a=commit;h=45abcec4b1ea764dec4bc56603feb9a951191388">AX_PROG_DOXYGEN: fix doxygen-ps rule</a></td>
<td class="link"><a href="/gitweb/?p=autoconf-archive.git;a=commit;h=45abcec4b1ea764dec4bc56603feb9a951191388">commit</a> | <a href="/gitweb/?p=autoconf-archive.git;a=commitdiff;h=45abcec4b1ea764dec4bc56603feb9a951191388">commitdiff</a> | <a href="/gitweb/?p=autoconf-archive.git;a=tree;h=45abcec4b1ea764dec4bc56603feb9a951191388;hb=45abcec4b1ea764dec4bc56603feb9a951191388">tree</a> | <a title="in format: tar.gz" href="/gitweb/?p=autoconf-archive.git;a=snapshot;h=45abcec4b1ea764dec4bc56603feb9a951191388;sf=tgz">snapshot</a></td>
</tr>
<tr class="light">
<td title="6 weeks ago"><i>2015-11-23</i></td>
<td class="author"><a title="Search for commits authored by Olaf Mandel" class="list" href="/gitweb/?p=autoconf-archive.git;a=search;s=Olaf+Mandel;st=author">Olaf Mandel</a></td><td><a class="list subject" href="/gitweb/?p=autoconf-archive.git;a=commit;h=9d9eb0f007c4ca2e1e0e9f0c6e3282a435a36fd3">AX_PROG_DOXYGEN: add support for silent rules</a></td>
<td class="link"><a href="/gitweb/?p=autoconf-archive.git;a=commit;h=9d9eb0f007c4ca2e1e0e9f0c6e3282a435a36fd3">commit</a> | <a href="/gitweb/?p=autoconf-archive.git;a=commitdiff;h=9d9eb0f007c4ca2e1e0e9f0c6e3282a435a36fd3">commitdiff</a> | <a href="/gitweb/?p=autoconf-archive.git;a=tree;h=9d9eb0f007c4ca2e1e0e9f0c6e3282a435a36fd3;hb=9d9eb0f007c4ca2e1e0e9f0c6e3282a435a36fd3">tree</a> | <a title="in format: tar.gz" href="/gitweb/?p=autoconf-archive.git;a=snapshot;h=9d9eb0f007c4ca2e1e0e9f0c6e3282a435a36fd3;sf=tgz">snapshot</a></td>
</tr>
<tr class="dark">
<td title="6 weeks ago"><i>2015-11-23</i></td>
<td class="author"><a title="Search for commits authored by Olaf Mandel" class="list" href="/gitweb/?p=autoconf-archive.git;a=search;s=Olaf+Mandel;st=author">Olaf Mandel</a></td><td><a class="list subject" href="/gitweb/?p=autoconf-archive.git;a=commit;h=5ef1e0b754692f089fa031100deeca1d05693fe6">AX_PROG_DOXYGEN: provide a DX_RULES substitution</a></td>
<td class="link"><a href="/gitweb/?p=autoconf-archive.git;a=commit;h=5ef1e0b754692f089fa031100deeca1d05693fe6">commit</a> | <a href="/gitweb/?p=autoconf-archive.git;a=commitdiff;h=5ef1e0b754692f089fa031100deeca1d05693fe6">commitdiff</a> | <a href="/gitweb/?p=autoconf-archive.git;a=tree;h=5ef1e0b754692f089fa031100deeca1d05693fe6;hb=5ef1e0b754692f089fa031100deeca1d05693fe6">tree</a> | <a title="in format: tar.gz" href="/gitweb/?p=autoconf-archive.git;a=snapshot;h=5ef1e0b754692f089fa031100deeca1d05693fe6;sf=tgz">snapshot</a></td>
</tr>
<tr class="light">
<td title="6 weeks ago"><i>2015-11-23</i></td>
<td class="author"><a title="Search for commits authored by Peter Simons" class="list" href="/gitweb/?p=autoconf-archive.git;a=search;s=Peter+Simons;st=author">Peter Simons</a></td><td><a class="list subject" href="/gitweb/?p=autoconf-archive.git;a=commit;h=0d7c7750e02dbd9ebe5373e470a14bcf29c1e2b9">Merge pull request #47 from emikulic/lcov</a></td>
<td class="link"><a href="/gitweb/?p=autoconf-archive.git;a=commit;h=0d7c7750e02dbd9ebe5373e470a14bcf29c1e2b9">commit</a> | <a href="/gitweb/?p=autoconf-archive.git;a=commitdiff;h=0d7c7750e02dbd9ebe5373e470a14bcf29c1e2b9">commitdiff</a> | <a href="/gitweb/?p=autoconf-archive.git;a=tree;h=0d7c7750e02dbd9ebe5373e470a14bcf29c1e2b9;hb=0d7c7750e02dbd9ebe5373e470a14bcf29c1e2b9">tree</a> | <a title="in format: tar.gz" href="/gitweb/?p=autoconf-archive.git;a=snapshot;h=0d7c7750e02dbd9ebe5373e470a14bcf29c1e2b9;sf=tgz">snapshot</a></td>
</tr>
<tr class="dark">
<td title="6 weeks ago"><i>2015-11-23</i></td>
<td class="author"><a title="Search for commits authored by Emil Mikulic" class="list" href="/gitweb/?p=autoconf-archive.git;a=search;s=Emil+Mikulic;st=author">Emil Mikulic</a></td><td><a class="list subject" href="/gitweb/?p=autoconf-archive.git;a=commit;h=b49708bdb5ba72e075c94b77daf0f361a2752e2e">Add 1.12 to list of lcov versions.</a></td>
<td class="link"><a href="/gitweb/?p=autoconf-archive.git;a=commit;h=b49708bdb5ba72e075c94b77daf0f361a2752e2e">commit</a> | <a href="/gitweb/?p=autoconf-archive.git;a=commitdiff;h=b49708bdb5ba72e075c94b77daf0f361a2752e2e">commitdiff</a> | <a href="/gitweb/?p=autoconf-archive.git;a=tree;h=b49708bdb5ba72e075c94b77daf0f361a2752e2e;hb=b49708bdb5ba72e075c94b77daf0f361a2752e2e">tree</a> | <a title="in format: tar.gz" href="/gitweb/?p=autoconf-archive.git;a=snapshot;h=b49708bdb5ba72e075c94b77daf0f361a2752e2e;sf=tgz">snapshot</a></td>
</tr>
<tr class="light">
<td title="2 months ago"><i>2015-10-13</i></td>
<td class="author"><a title="Search for commits authored by Olaf Mandel" class="list" href="/gitweb/?p=autoconf-archive.git;a=search;s=Olaf+Mandel;st=author">Olaf Mandel</a></td><td><a class="list subject" href="/gitweb/?p=autoconf-archive.git;a=commit;h=7fadb5b146f81e90d661dc89ad0edc1e4d239dd6">AX_CODE_COVERAGE: add CPPFLAGS to skip assertions</a></td>
<td class="link"><a href="/gitweb/?p=autoconf-archive.git;a=commit;h=7fadb5b146f81e90d661dc89ad0edc1e4d239dd6">commit</a> | <a href="/gitweb/?p=autoconf-archive.git;a=commitdiff;h=7fadb5b146f81e90d661dc89ad0edc1e4d239dd6">commitdiff</a> | <a href="/gitweb/?p=autoconf-archive.git;a=tree;h=7fadb5b146f81e90d661dc89ad0edc1e4d239dd6;hb=7fadb5b146f81e90d661dc89ad0edc1e4d239dd6">tree</a> | <a title="in format: tar.gz" href="/gitweb/?p=autoconf-archive.git;a=snapshot;h=7fadb5b146f81e90d661dc89ad0edc1e4d239dd6;sf=tgz">snapshot</a></td>
</tr>
<tr>
<td colspan="4"><a href="/gitweb/?p=autoconf-archive.git;a=shortlog">...</a></td>
</tr>
</table>
<div class="header">
<a class="title" href="/gitweb/?p=autoconf-archive.git;a=tags">tags</a>
</div>
<table class="tags">
<tr class="dark">
<td><i>3 months ago</i></td>
<td><a class="list name" href="/gitweb/?p=autoconf-archive.git;a=commit;h=557b34d93d5a03fd66a2d40b02152b8f12eccf4c">v2015.09.25</a></td>
<td><a class="list subject" title="GNU Autoconf Archive Version 2015.09.25" href="/gitweb/?p=autoconf-archive.git;a=tag;h=501f40aa459c64f2cb2ccf22b16e56d37d718b57">GNU Autoconf Archive Version 2015... </a></td>
<td class="selflink"><a href="/gitweb/?p=autoconf-archive.git;a=tag;h=501f40aa459c64f2cb2ccf22b16e56d37d718b57">tag</a></td>
<td class="link"> | <a href="/gitweb/?p=autoconf-archive.git;a=commit;h=557b34d93d5a03fd66a2d40b02152b8f12eccf4c">commit</a> | <a href="/gitweb/?p=autoconf-archive.git;a=shortlog;h=refs/tags/v2015.09.25">shortlog</a> | <a href="/gitweb/?p=autoconf-archive.git;a=log;h=refs/tags/v2015.09.25">log</a></td>
</tr><tr class="light">
<td><i>10 months ago</i></td>
<td><a class="list name" href="/gitweb/?p=autoconf-archive.git;a=commit;h=8c3a792472f405e3b89f9a2a48a700caa67089d8">v2015.02.24</a></td>
<td><a class="list subject" href="/gitweb/?p=autoconf-archive.git;a=tag;h=927daaca51b7c8a67f5a1f0fbca104b8316e92fd">Autoconf Archive Version 2015.02.24</a></td>
<td class="selflink"><a href="/gitweb/?p=autoconf-archive.git;a=tag;h=927daaca51b7c8a67f5a1f0fbca104b8316e92fd">tag</a></td>
<td class="link"> | <a href="/gitweb/?p=autoconf-archive.git;a=commit;h=8c3a792472f405e3b89f9a2a48a700caa67089d8">commit</a> | <a href="/gitweb/?p=autoconf-archive.git;a=shortlog;h=refs/tags/v2015.02.24">shortlog</a> | <a href="/gitweb/?p=autoconf-archive.git;a=log;h=refs/tags/v2015.02.24">log</a></td>
</tr><tr class="dark">
<td><i>11 months ago</i></td>
<td><a class="list name" href="/gitweb/?p=autoconf-archive.git;a=commit;h=6ddf085ffd6c699070fe61cdc09e76743a49875f">v2015.02.04</a></td>
<td><a class="list subject" href="/gitweb/?p=autoconf-archive.git;a=tag;h=06ab78be731e434efe39c11a91a70df7b527c386">Autoconf Archive Version 2015.02.04</a></td>
<td class="selflink"><a href="/gitweb/?p=autoconf-archive.git;a=tag;h=06ab78be731e434efe39c11a91a70df7b527c386">tag</a></td>
<td class="link"> | <a href="/gitweb/?p=autoconf-archive.git;a=commit;h=6ddf085ffd6c699070fe61cdc09e76743a49875f">commit</a> | <a href="/gitweb/?p=autoconf-archive.git;a=shortlog;h=refs/tags/v2015.02.04">shortlog</a> | <a href="/gitweb/?p=autoconf-archive.git;a=log;h=refs/tags/v2015.02.04">log</a></td>
</tr><tr class="light">
<td><i>14 months ago</i></td>
<td><a class="list name" href="/gitweb/?p=autoconf-archive.git;a=commit;h=7f56aaabdac6dc45ee083c6f3afdfc412a600124">v2014.10.15</a></td>
<td><a class="list subject" href="/gitweb/?p=autoconf-archive.git;a=tag;h=34f91518355f1f2cba082678c60008a1ed8d009a">Autoconf Archive Version 2014.10.15</a></td>
<td class="selflink"><a href="/gitweb/?p=autoconf-archive.git;a=tag;h=34f91518355f1f2cba082678c60008a1ed8d009a">tag</a></td>
<td class="link"> | <a href="/gitweb/?p=autoconf-archive.git;a=commit;h=7f56aaabdac6dc45ee083c6f3afdfc412a600124">commit</a> | <a href="/gitweb/?p=autoconf-archive.git;a=shortlog;h=refs/tags/v2014.10.15">shortlog</a> | <a href="/gitweb/?p=autoconf-archive.git;a=log;h=refs/tags/v2014.10.15">log</a></td>
</tr><tr class="dark">
<td><i>22 months ago</i></td>
<td><a class="list name" href="/gitweb/?p=autoconf-archive.git;a=commit;h=8935e464cac3a8eb584c093b35863ba08690471b">v2014.02.28</a></td>
<td><a class="list subject" title="GNU Autoconf Archive Version 2014.02.28" href="/gitweb/?p=autoconf-archive.git;a=tag;h=42b15b27650f65c19a707bee7db2f987ce69a0dd">GNU Autoconf Archive Version 2014... </a></td>
<td class="selflink"><a href="/gitweb/?p=autoconf-archive.git;a=tag;h=42b15b27650f65c19a707bee7db2f987ce69a0dd">tag</a></td>
<td class="link"> | <a href="/gitweb/?p=autoconf-archive.git;a=commit;h=8935e464cac3a8eb584c093b35863ba08690471b">commit</a> | <a href="/gitweb/?p=autoconf-archive.git;a=shortlog;h=refs/tags/v2014.02.28">shortlog</a> | <a href="/gitweb/?p=autoconf-archive.git;a=log;h=refs/tags/v2014.02.28">log</a></td>
</tr><tr class="light">
<td><i>2 years ago</i></td>
<td><a class="list name" href="/gitweb/?p=autoconf-archive.git;a=commit;h=dca0d5c7de9729a45ebb8934db114912dc04e37b">v2013.11.01</a></td>
<td><a class="list subject" title="GNU Autoconf Archive Version 2013.11.01" href="/gitweb/?p=autoconf-archive.git;a=tag;h=7c9026664523434623f3cc26d9e82dee00686e86">GNU Autoconf Archive Version 2013... </a></td>
<td class="selflink"><a href="/gitweb/?p=autoconf-archive.git;a=tag;h=7c9026664523434623f3cc26d9e82dee00686e86">tag</a></td>
<td class="link"> | <a href="/gitweb/?p=autoconf-archive.git;a=commit;h=dca0d5c7de9729a45ebb8934db114912dc04e37b">commit</a> | <a href="/gitweb/?p=autoconf-archive.git;a=shortlog;h=refs/tags/v2013.11.01">shortlog</a> | <a href="/gitweb/?p=autoconf-archive.git;a=log;h=refs/tags/v2013.11.01">log</a></td>
</tr><tr class="dark">
<td><i>2 years ago</i></td>
<td><a class="list name" href="/gitweb/?p=autoconf-archive.git;a=commit;h=f5f6d113d94f121210370e10f41045491581bdc2">v2013.06.09</a></td>
<td><a class="list subject" title="GNU Autoconf Archive Version 2013.06.09" href="/gitweb/?p=autoconf-archive.git;a=tag;h=7167a25793d36ad25d7068c2625d4cad78a83b26">GNU Autoconf Archive Version 2013... </a></td>
<td class="selflink"><a href="/gitweb/?p=autoconf-archive.git;a=tag;h=7167a25793d36ad25d7068c2625d4cad78a83b26">tag</a></td>
<td class="link"> | <a href="/gitweb/?p=autoconf-archive.git;a=commit;h=f5f6d113d94f121210370e10f41045491581bdc2">commit</a> | <a href="/gitweb/?p=autoconf-archive.git;a=shortlog;h=refs/tags/v2013.06.09">shortlog</a> | <a href="/gitweb/?p=autoconf-archive.git;a=log;h=refs/tags/v2013.06.09">log</a></td>
</tr><tr class="light">
<td><i>2 years ago</i></td>
<td><a class="list name" href="/gitweb/?p=autoconf-archive.git;a=commit;h=60e854d960fb96105bbe59449b21aa6a9d196de3">v2013.04.06</a></td>
<td><a class="list subject" title="GNU Autoconf Archive Version 2013.04.06" href="/gitweb/?p=autoconf-archive.git;a=tag;h=1d216b9587c50f9d6e59ad3bec8f8b0efd9242c5">GNU Autoconf Archive Version 2013... </a></td>
<td class="selflink"><a href="/gitweb/?p=autoconf-archive.git;a=tag;h=1d216b9587c50f9d6e59ad3bec8f8b0efd9242c5">tag</a></td>
<td class="link"> | <a href="/gitweb/?p=autoconf-archive.git;a=commit;h=60e854d960fb96105bbe59449b21aa6a9d196de3">commit</a> | <a href="/gitweb/?p=autoconf-archive.git;a=shortlog;h=refs/tags/v2013.04.06">shortlog</a> | <a href="/gitweb/?p=autoconf-archive.git;a=log;h=refs/tags/v2013.04.06">log</a></td>
</tr><tr class="dark">
<td><i>2 years ago</i></td>
<td><a class="list name" href="/gitweb/?p=autoconf-archive.git;a=commit;h=39acabca5d882efa48dcc6954af87e34bc11f563">v2013.02.02</a></td>
<td><a class="list subject" title="GNU Autoconf Archive Version 2013.02.02" href="/gitweb/?p=autoconf-archive.git;a=tag;h=8a60961910acc03aef17518b92924aeadeee4f51">GNU Autoconf Archive Version 2013... </a></td>
<td class="selflink"><a href="/gitweb/?p=autoconf-archive.git;a=tag;h=8a60961910acc03aef17518b92924aeadeee4f51">tag</a></td>
<td class="link"> | <a href="/gitweb/?p=autoconf-archive.git;a=commit;h=39acabca5d882efa48dcc6954af87e34bc11f563">commit</a> | <a href="/gitweb/?p=autoconf-archive.git;a=shortlog;h=refs/tags/v2013.02.02">shortlog</a> | <a href="/gitweb/?p=autoconf-archive.git;a=log;h=refs/tags/v2013.02.02">log</a></td>
</tr><tr class="light">
<td><i>3 years ago</i></td>
<td><a class="list name" href="/gitweb/?p=autoconf-archive.git;a=commit;h=46143062e93089a64550493e3a9659117559c662">v2012.11.14</a></td>
<td><a class="list subject" title="GNU Autoconf Archive Version 2012.11.14" href="/gitweb/?p=autoconf-archive.git;a=tag;h=627611bcad3d47691076ce5240760d44d41d95d2">GNU Autoconf Archive Version 2012... </a></td>
<td class="selflink"><a href="/gitweb/?p=autoconf-archive.git;a=tag;h=627611bcad3d47691076ce5240760d44d41d95d2">tag</a></td>
<td class="link"> | <a href="/gitweb/?p=autoconf-archive.git;a=commit;h=46143062e93089a64550493e3a9659117559c662">commit</a> | <a href="/gitweb/?p=autoconf-archive.git;a=shortlog;h=refs/tags/v2012.11.14">shortlog</a> | <a href="/gitweb/?p=autoconf-archive.git;a=log;h=refs/tags/v2012.11.14">log</a></td>
</tr><tr class="dark">
<td><i>3 years ago</i></td>
<td><a class="list name" href="/gitweb/?p=autoconf-archive.git;a=commit;h=6f4975c5cbfcf9c71852e75596b4b74f88c83195">v2012.09.08</a></td>
<td><a class="list subject" title="GNU Autoconf Archive Version 2012.09.08" href="/gitweb/?p=autoconf-archive.git;a=tag;h=42f4f4d30a2c03df0ab991a4391cb336a1aa408f">GNU Autoconf Archive Version 2012... </a></td>
<td class="selflink"><a href="/gitweb/?p=autoconf-archive.git;a=tag;h=42f4f4d30a2c03df0ab991a4391cb336a1aa408f">tag</a></td>
<td class="link"> | <a href="/gitweb/?p=autoconf-archive.git;a=commit;h=6f4975c5cbfcf9c71852e75596b4b74f88c83195">commit</a> | <a href="/gitweb/?p=autoconf-archive.git;a=shortlog;h=refs/tags/v2012.09.08">shortlog</a> | <a href="/gitweb/?p=autoconf-archive.git;a=log;h=refs/tags/v2012.09.08">log</a></td>
</tr><tr class="light">
<td><i>3 years ago</i></td>
<td><a class="list name" href="/gitweb/?p=autoconf-archive.git;a=commit;h=d94f9f602747b87f97912fa649144bed4339839d">v2012.04.07</a></td>
<td><a class="list subject" title="GNU Autoconf Archive Version 2012.04.07" href="/gitweb/?p=autoconf-archive.git;a=tag;h=a2b0e3d1b50469dcc98d4c7f23933fb549effdac">GNU Autoconf Archive Version 2012... </a></td>
<td class="selflink"><a href="/gitweb/?p=autoconf-archive.git;a=tag;h=a2b0e3d1b50469dcc98d4c7f23933fb549effdac">tag</a></td>
<td class="link"> | <a href="/gitweb/?p=autoconf-archive.git;a=commit;h=d94f9f602747b87f97912fa649144bed4339839d">commit</a> | <a href="/gitweb/?p=autoconf-archive.git;a=shortlog;h=refs/tags/v2012.04.07">shortlog</a> | <a href="/gitweb/?p=autoconf-archive.git;a=log;h=refs/tags/v2012.04.07">log</a></td>
</tr><tr class="dark">
<td><i>4 years ago</i></td>
<td><a class="list name" href="/gitweb/?p=autoconf-archive.git;a=commit;h=77bb2ab18a583eca3717406e47cba73426f59ab1">v2011.12.21</a></td>
<td><a class="list subject" title="GNU Autoconf Archive Version 2011.12.21" href="/gitweb/?p=autoconf-archive.git;a=tag;h=b57da13e32c04d5a800257ab1b3ca48cd7aa1871">GNU Autoconf Archive Version 2011... </a></td>
<td class="selflink"><a href="/gitweb/?p=autoconf-archive.git;a=tag;h=b57da13e32c04d5a800257ab1b3ca48cd7aa1871">tag</a></td>
<td class="link"> | <a href="/gitweb/?p=autoconf-archive.git;a=commit;h=77bb2ab18a583eca3717406e47cba73426f59ab1">commit</a> | <a href="/gitweb/?p=autoconf-archive.git;a=shortlog;h=refs/tags/v2011.12.21">shortlog</a> | <a href="/gitweb/?p=autoconf-archive.git;a=log;h=refs/tags/v2011.12.21">log</a></td>
</tr><tr class="light">
<td><i>4 years ago</i></td>
<td><a class="list name" href="/gitweb/?p=autoconf-archive.git;a=commit;h=e6d6aa7c38f235c09efbe97976454e4f65a52539">v2011.09.17</a></td>
<td><a class="list subject" title="GNU Autoconf Archive Version 2011.09.17" href="/gitweb/?p=autoconf-archive.git;a=tag;h=7d6e0f275e488dc135ac0b73330ec69978f9dabf">GNU Autoconf Archive Version 2011... </a></td>
<td class="selflink"><a href="/gitweb/?p=autoconf-archive.git;a=tag;h=7d6e0f275e488dc135ac0b73330ec69978f9dabf">tag</a></td>
<td class="link"> | <a href="/gitweb/?p=autoconf-archive.git;a=commit;h=e6d6aa7c38f235c09efbe97976454e4f65a52539">commit</a> | <a href="/gitweb/?p=autoconf-archive.git;a=shortlog;h=refs/tags/v2011.09.17">shortlog</a> | <a href="/gitweb/?p=autoconf-archive.git;a=log;h=refs/tags/v2011.09.17">log</a></td>
</tr><tr class="dark">
<td><i>4 years ago</i></td>
<td><a class="list name" href="/gitweb/?p=autoconf-archive.git;a=commit;h=756c399edef5a71f73f05b61cefde0ba4842e6c5">v2011.07.17</a></td>
<td><a class="list subject" title="GNU Autoconf Archive Version 2011.07.17" href="/gitweb/?p=autoconf-archive.git;a=tag;h=7630d0be0daee35103e690dcd98a41e7780dbad9">GNU Autoconf Archive Version 2011... </a></td>
<td class="selflink"><a href="/gitweb/?p=autoconf-archive.git;a=tag;h=7630d0be0daee35103e690dcd98a41e7780dbad9">tag</a></td>
<td class="link"> | <a href="/gitweb/?p=autoconf-archive.git;a=commit;h=756c399edef5a71f73f05b61cefde0ba4842e6c5">commit</a> | <a href="/gitweb/?p=autoconf-archive.git;a=shortlog;h=refs/tags/v2011.07.17">shortlog</a> | <a href="/gitweb/?p=autoconf-archive.git;a=log;h=refs/tags/v2011.07.17">log</a></td>
</tr><tr class="light">
<td><i>4 years ago</i></td>
<td><a class="list name" href="/gitweb/?p=autoconf-archive.git;a=commit;h=7fd27beef29f4d6ff69f4fce874185522174806b">v2011.04.12</a></td>
<td><a class="list subject" title="GNU Autoconf Archive Version 2011.04.12" href="/gitweb/?p=autoconf-archive.git;a=tag;h=899801493b5120e9d1aaa75146adf04acff3cc99">GNU Autoconf Archive Version 2011... </a></td>
<td class="selflink"><a href="/gitweb/?p=autoconf-archive.git;a=tag;h=899801493b5120e9d1aaa75146adf04acff3cc99">tag</a></td>
<td class="link"> | <a href="/gitweb/?p=autoconf-archive.git;a=commit;h=7fd27beef29f4d6ff69f4fce874185522174806b">commit</a> | <a href="/gitweb/?p=autoconf-archive.git;a=shortlog;h=refs/tags/v2011.04.12">shortlog</a> | <a href="/gitweb/?p=autoconf-archive.git;a=log;h=refs/tags/v2011.04.12">log</a></td>
</tr><tr>
<td colspan="5"><a href="/gitweb/?p=autoconf-archive.git;a=tags">...</a></td>
</tr>
</table>
<div class="header">
<a class="title" href="/gitweb/?p=autoconf-archive.git;a=heads">heads</a>
</div>
<table class="heads">
<tr class="dark">
<td><i>3 weeks ago</i></td>
<td class="current_head"><a class="list name" href="/gitweb/?p=autoconf-archive.git;a=shortlog;h=refs/heads/master">master</a></td>
<td class="link"><a href="/gitweb/?p=autoconf-archive.git;a=shortlog;h=refs/heads/master">shortlog</a> | <a href="/gitweb/?p=autoconf-archive.git;a=log;h=refs/heads/master">log</a> | <a href="/gitweb/?p=autoconf-archive.git;a=tree;h=refs/heads/master;hb=master">tree</a></td>
</tr><tr class="light">
<td><i>6 years ago</i></td>
<td><a class="list name" href="/gitweb/?p=autoconf-archive.git;a=shortlog;h=refs/heads/pre-savannah-history">pre-savannah-history</a></td>
<td class="link"><a href="/gitweb/?p=autoconf-archive.git;a=shortlog;h=refs/heads/pre-savannah-history">shortlog</a> | <a href="/gitweb/?p=autoconf-archive.git;a=log;h=refs/heads/pre-savannah-history">log</a> | <a href="/gitweb/?p=autoconf-archive.git;a=tree;h=refs/heads/pre-savannah-history;hb=pre-savannah-history">tree</a></td>
</tr></table>
<div class="page_footer">
<div class="page_footer_text">GNU Autoconf Archive</div>
<a class="rss_logo" title="log RSS feed" href="/gitweb/?p=autoconf-archive.git;a=rss">RSS</a>
<a class="rss_logo" title="log Atom feed" href="/gitweb/?p=autoconf-archive.git;a=atom">Atom</a>
</div>
<script type="text/javascript" src="gitweb.js"></script>
</body>
</html>