/**
 * export_pdf.mjs — Windows Chrome CDP 経由で slides.html → slides.pdf を生成する
 * 実行: node export_pdf.mjs
 *
 * 仕組み:
 *  1. Windows Chrome を --remote-debugging-address=0.0.0.0 で起動
 *  2. WSL から Windows ホスト IP 経由で Chrome DevTools Protocol に接続
 *  3. ?print-pdf URL（白背景モード）を開き、Reveal.js 初期化待機後に PDF 出力
 *  4. base64 で受け取った PDF データをファイルに書き込む
 */

import { spawn } from 'child_process';
import http from 'http';
import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';
import { createRequire } from 'module';

const require = createRequire(import.meta.url);
const WebSocket = require('ws');

const __dirname  = path.dirname(fileURLToPath(import.meta.url));

// ── パス設定 ──────────────────────────────────────────────────────
const CHROME      = '/mnt/c/Program Files/Google/Chrome/Application/chrome.exe';
const DEBUG_PORT  = 9221;

// WSL 側 Linux パス → Windows ファイル URL へ変換
// /mnt/d/foo/bar → file:///D:/foo/bar
function toFileUrl(linuxPath) {
  return 'file:///' + linuxPath
    .replace(/^\/mnt\/([a-z])\//, (_, d) => d.toUpperCase() + ':/')
    .replace(/\//g, '/');
}

const SLIDE_LINUX  = path.join(__dirname, 'slides.html');
const OUTPUT_LINUX = path.join(__dirname, 'slides.pdf');
const SLIDE_URL    = toFileUrl(SLIDE_LINUX) + '?print-pdf';

console.log('Slide URL:', SLIDE_URL);
console.log('Output   :', OUTPUT_LINUX);

// ── Windows ホスト IP 取得 ───────────────────────────────────────
function getWindowsHostIp() {
  try {
    const conf = fs.readFileSync('/etc/resolv.conf', 'utf8');
    const m = conf.match(/nameserver\s+(\S+)/);
    return m ? m[1] : '172.31.224.1';
  } catch { return '172.31.224.1'; }
}
const WIN_HOST = getWindowsHostIp();
console.log('Windows host IP:', WIN_HOST);

// ── Chrome 起動 ──────────────────────────────────────────────────
let chromeProcess = null;

function launchChrome() {
  const userDataDir = `/tmp/chrome-pdf-${Date.now()}`;
  chromeProcess = spawn(CHROME, [
    `--remote-debugging-port=${DEBUG_PORT}`,
    '--remote-debugging-address=0.0.0.0',   // WSL からアクセスできるよう全 IF でリスン
    '--headless=new',
    '--no-sandbox',
    '--disable-gpu',
    '--disable-setuid-sandbox',
    '--no-first-run',
    '--no-default-browser-check',
    '--allow-file-access-from-files',        // file:// から相対パス PNG を読む
    '--disable-web-security',               // file:// の CORS を無効
    `--user-data-dir=${userDataDir}`,
  ], { stdio: 'ignore' });

  chromeProcess.on('error', (e) => { console.error('Chrome error:', e.message); });
  return chromeProcess;
}

// ── CDP JSON 取得 ─────────────────────────────────────────────────
function getCdpVersion(host, port, retries = 20) {
  return new Promise((resolve, reject) => {
    const attempt = (n) => {
      http.get(`http://${host}:${port}/json/version`, (res) => {
        let d = '';
        res.on('data', c => d += c);
        res.on('end', () => {
          try { resolve(JSON.parse(d)); }
          catch { if (n > 0) setTimeout(() => attempt(n - 1), 800); else reject(new Error('CDP parse error')); }
        });
      }).on('error', () => {
        if (n > 0) setTimeout(() => attempt(n - 1), 800);
        else reject(new Error(`Chrome CDP not responding at ${host}:${port}`));
      });
    };
    attempt(retries);
  });
}

// ── CDP セッション ────────────────────────────────────────────────
function cdp(ws) {
  let id = 1;
  const pending = new Map();
  ws.on('message', (raw) => {
    const msg = JSON.parse(raw);
    if (msg.id && pending.has(msg.id)) {
      const { res, rej } = pending.get(msg.id);
      pending.delete(msg.id);
      msg.error ? rej(new Error(msg.error.message)) : res(msg.result);
    }
  });
  return {
    send(method, params = {}) {
      return new Promise((res, rej) => {
        const sid = id++;
        pending.set(sid, { res, rej });
        ws.send(JSON.stringify({ id: sid, method, params }));
      });
    },
    once(event) {
      return new Promise((resolve) => {
        const h = (raw) => {
          const msg = JSON.parse(raw);
          if (msg.method === event) { ws.off('message', h); resolve(msg.params); }
        };
        ws.on('message', h);
      });
    }
  };
}

const delay = (ms) => new Promise(r => setTimeout(r, ms));

// ── メイン ────────────────────────────────────────────────────────
(async () => {
  launchChrome();
  console.log('Chrome launched, waiting for CDP...');

  let version;
  try {
    version = await getCdpVersion(WIN_HOST, DEBUG_PORT);
  } catch (e) {
    console.error('Failed to connect Chrome CDP:', e.message);
    chromeProcess?.kill(); process.exit(1);
  }

  // WSL → Windows IP に変換した WS URL
  const wsUrl = version.webSocketDebuggerUrl.replace(/127\.0\.0\.1/, WIN_HOST);
  console.log('WS:', wsUrl);

  const ws = new WebSocket(wsUrl);
  ws.on('error', (e) => { console.error('WS error:', e.message); chromeProcess?.kill(); process.exit(1); });

  await new Promise(r => ws.once('open', r));
  const session = cdp(ws);

  // ページナビゲーション
  await session.send('Page.enable');
  await session.send('Runtime.enable');

  const loaded = session.once('Page.loadEventFired');
  await session.send('Page.navigate', { url: SLIDE_URL });
  console.log('Navigating...');

  await Promise.race([loaded, delay(30000)]);
  console.log('Page loaded. Waiting for Reveal.js + MathJax...');

  // Reveal.js 準備完了を待つ
  await session.send('Runtime.evaluate', {
    expression: `new Promise(resolve => {
      const t = setInterval(() => {
        if (window.Reveal && Reveal.isReady && Reveal.isReady()) {
          clearInterval(t); resolve(true);
        }
      }, 500);
      setTimeout(() => { clearInterval(t); resolve(false); }, 20000);
    })`,
    awaitPromise: true,
    timeout: 25000,
  });
  console.log('Reveal.js ready.');

  // MathJax・フォント追加待機
  await delay(5000);
  console.log('Printing to PDF...');

  // 1280×720 px → inch 換算（96dpi 基準）
  const { data } = await session.send('Page.printToPDF', {
    printBackground : true,
    paperWidth      : 1280 / 96,  // ≈ 13.33 inch
    paperHeight     : 720  / 96,  // = 7.5 inch
    marginTop       : 0,
    marginBottom    : 0,
    marginLeft      : 0,
    marginRight     : 0,
    scale           : 1.0,
  });

  const buf = Buffer.from(data, 'base64');
  fs.writeFileSync(OUTPUT_LINUX, buf);
  console.log(`✅ Done: ${OUTPUT_LINUX} (${(buf.length / 1024).toFixed(0)} KB)`);

  ws.close();
  chromeProcess?.kill();
  process.exit(0);
})();
